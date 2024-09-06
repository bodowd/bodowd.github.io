# duckdb_rdkit

Experiments in building a cheminformatics extension for an OLAP database

---

# Background information

- Two main types of database workloads
  <br>

  - Online Transactional Processing (OLTP) & Online Analytical Processing (OLAP)
    <br>
    - access few or single row(s) vs access many or all rows
      <br>
      - update amount in an account vs average account balance
        <br>
    - the techniques used in building these systems are very different:
      - OLTP optimizes for finding individual rows
        <br>
      - OLAP optimizes for full column scans
        <br>
      - i.e. different strategies to indexing, storage formats, query execution

- duckdb is an in-process OLAP database and query execution engine
  <br>
- RDKit is a cheminformatics toolkit for working with molecules on computers

---

# Implementing the duckdb_rdkit extension

## Exact match and substructure match

- is this molecule the same as that molecule?
  <br>
- is this fragment found in that molecule?

## The Mol object

```
~~~graph-easy --as=boxart
[ CCO (SMILES) ]  <-> [ RDKit Mol object (in-memory) ]  <-> [ Serialize to binary (on-disk)]
~~~
```

- The RDKit Mol object contains a lot of information about the molecule for
  cheminformatics work

## Initial attempt at exact match in duckdb_rdkit:

- adapted the code for comparing molecules from RDKit extensions for Postgres & SQLite
  <br>
- very poor performance

| Query | Standard method (s) |
| :---- | :------------------ |
| 1     | 17.238              |
| 2     | 12.555              |
| 3     | 22.196              |
| 4     | 12.245              |

- chembl 33, ~2.3 million molecules
- default duckdb settings: AMD Ryzen 5 4500U CPU, 16GB RAM, Samsung PM991 SSD

---

# Implementing the duckdb_rdkit extension

## Improving exact match and substructure match performance

- Fast forward after many experiments and iterations...
  <br>

- Inspired by Umbra-style strings (Neumann, T., Freitag M. CIDR 2020)
  <br>
- "Umbra Mol"

---

# Implementing the duckdb_rdkit extension

## Short-circuiting for exact match

- Examined the code for comparing molecules in the Postgres and SQLite RDKit extensions

```c++
bool mol_cmp(const RDKit::ROMol &m1, const RDKit::ROMol &m2) {
  int res = m1.getNumAtoms() - m2.getNumAtoms();

  if (res) {
    return false;
  }

  res = m1.getNumBonds() - m2.getNumBonds();
  if (res) {
    return false;
  }

  res = int(RDKit::Descriptors::calcAMW(m1, false) -
            RDKit::Descriptors::calcAMW(m2, false) + .5);
  if (res) {
    return false;
  }

  res = m1.getRingInfo()->numRings() - m2.getRingInfo()->numRings();
  if (res) {
    return false;
  }

  // more expensive checks like substructure match and creating canonical
  // smiles. Omitted for brevity

}
```

- short-circuits with cheaper comparisons
  <br>
- but requires deserialization of the binary to the in-memory RDKit Mol

---

# Implementing the duckdb_rdkit extension

## Short-circuiting for exact match

- Store the counts in a prefix in front of the binary molecule

```
~~~graph-easy --as=boxart
[num_atoms | num_bonds | amw | num_rings | binary rdkit molecule]
~~~
```

- Compare prefix binary directly
  <br>
  - return false if prefixes don't match
    <br>
  - no deserialization to RDKit Mol
    <br>
  - if prefixes match, deserialize and do the rest of the checks

```
~~~graph-easy --as=boxart
[ CCO (SMILES) ]  <-> [ RDKit Mol object (in-memory) ]  <-> [ Serialize to binary (on-disk)]
~~~
```

---

# Implementing the duckdb_rdkit extension

## Short-circuiting for exact match

- Analyzed chembl 33 to optimize size of prefix
  <br>
  - in some cases, 99% of the data can be represented with a few bits (initial experiments I used 4B each)

```
D select quantile_cont(num_atoms, 0.99), quantile_cont(num_bonds, 0.99), quantile_cont(num_rings, 0.99), quantile_cont(amw::integer, 0.99) from molecule;
┌────────────────────────────────┬────────────────────────────────┬────────────────────────────────┬───────────────────────────────────────────┐
│ quantile_cont(num_atoms, 0.99) │ quantile_cont(num_bonds, 0.99) │ quantile_cont(num_rings, 0.99) │ quantile_cont(CAST(amw AS INTEGER), 0.99) │
│             double             │             double             │             double             │                  double                   │
├────────────────────────────────┼────────────────────────────────┼────────────────────────────────┼───────────────────────────────────────────┤
│                          200.0 │                           37.0 │                            8.0 │                                    1446.0 │
└────────────────────────────────┴────────────────────────────────┴────────────────────────────────┴───────────────────────────────────────────┘
Run Time (s): real 0.190 user 0.182839 sys 0.164021

```

- number of atoms, 7 bits (oops - but still good: 97th percentile is 119 atoms)
- number of bonds, 6 bits
- number of rings, 3 bits
- amw, 11 bits

total: 27 bits (~4B as opposed to 20B)

---

# Implementing the duckdb_rdkit extension

## Short-circuiting for exact match

| Query | Standard method (s) | Umbra-mol (s) | speedup |
| :---- | :------------------ | :------------ | ------: |
| 1     | 17.238              | 0.496         |  34.75x |
| 2     | 12.555              | 0.473         |  26.54x |
| 3     | 22.196              | 0.364         |  60.98x |
| 4     | 12.245              | 0.359         |  34.11x |

- Avoiding deserialization is very good

---

# Implementing the duckdb_rdkit extension

## Substructure filter: dalke_fp

- Substructure comparison is computationally expensive: complicated graph stuff
  <br>
- Substructure filter:
  <br>

  - Query: NC, Target: OC -- not a substructure
    <br>
  - Query: CC, Target: NCC -- could be substructure, need isomorphism test
    <br>
  - Can only bail out for the first case: If the query has the feature, but the
    target does not, cannot be a substructure match
    <br>
  - Need a good set of features

- Andrew Dalke shared his work on developing a substructure filter (2012):
  <br>

  - Long story short: 55 bit fingerprint (`dalke_fp`) (1 if feature present, 0 if not) that
    we can use to short-circuit by checking if the bit is on in the query, but
    off in the target
    <br>

- Incorporated into the prefix of Umbra Mol for speeding up substructure matches
  <br>
- Details here:
  - http://www.dalkescientific.com/writings/diary/archive/2012/06/11/optimizing_substructure_keys.html
  - https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02078.html

---

# Implementing the duckdb_rdkit extension

## Store pointer to binary molecule instead of inlining

- During experiments with dalke_fp substructure filter, saw performance regression
  <br>
  - Tracked it down to the struct getting too big

```
~~~graph-easy --as=boxart
[counts (4B) | dalke_fp (8B)| binary rdkit molecule (chembl 33 avg: ~455B)]
~~~
```

- Storing a pointer to the binary molecule can shrink the struct significantly

```
~~~graph-easy --as=boxart
[counts (4B) | dalke_fp (8B)| pointer to binary molecule (8B -- 64 bit CPU)]
~~~
```

- Pointer swizzling is required (virtual memory address vs disk offset)
  <br>
  - Long story short: Made it look like a `string_t` in duckdb and duckdb's internals handles the rest

---

# Implementing the duckdb_rdkit extension

## Umbra-Mol results

```
~~~graph-easy --as=boxart
[counts (4B) | dalke_fp (8B)| pointer to binary molecule (8B -- 64 bit CPU)]
~~~
```

- Default duckdb & Postgres settings
- AMD Ryzen 5 4500U CPU, 16GB RAM, Samsung PM991 SSD
- Postgres running in docker
  - gist index on molecules in Postgres

##### Exact match

| Query | Standard method (s) | Umbra-mol v2 (s) | speedup (Umbra-mol vs standard method) | Postgres control (s) |
| :---- | :------------------ | :--------------- | :------------------------------------- | :------------------- |
| 1     | 17.238              | 0.179            | 96x                                    | 0.084                |
| 2     | 12.555              | 0.145            | 87x                                    | 233                  |
| 3     | 13.027              | 0.263            | 50x                                    | 2.47                 |
| 4     | 12.245              | 0.255            | 48x                                    | 6.185                |

##### Substructure match

| Query | Standard method (s) | Umbra-mol v2 (s) | speedup (Umbra-mol vs standard method) | Postgres control (s) |
| :---- | :------------------ | :--------------- | :------------------------------------- | :------------------- |
| 1     | 23.388              | 0.267            | 88x                                    | 0.741                |
| 2     | 14.094              | 5.93             | 2x                                     | 98                   |
| 3     | 14.294              | 0.553            | 26x                                    | 12.114               |
| 4     | 13.994              | 6.804            | 2x                                     | 1237 (20 min)        |

---

# Wrapping up

- Currently, not many features, and very experimental
- github.com/bodowd/duckdb_rdkit

---

# Acknowledgments

- Thomas Leyer
- Fuad Abdallah
- Andrew Dalke

bing.odowd@bayer.com
