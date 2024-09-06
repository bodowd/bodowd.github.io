# duckdb_rdkit: RDKit extension in duckdb

- Background:

  - OLTP vs OLAP
  - What is duckdb, and what is interesting about it?
  - Cheminformatics workloads & RDKit

- Implementation of duckdb_rdkit:
  - implementing the Mol type: some kind of format (i.e. SMILES) -> RDKit Mol Object (in memory representation) -> Serialize to binary (disk format)
  - implementing exact match
    - initial attempt
    - umbra_mol initial attempt (prefix, but inlined)
      - introduce umbra style string idea
      - results of intial attempt and the prefix
    - umbra_mol second attempt (prefix, pointer)
      - chembl analysis to improve the count prefix for exact match
      - dalke_fp for substructure match
      - pointer to the binary molecule
        - analyzed duckdb, where does it do pointer swizzling? string_t, etc.

---

# Very brief intro to OLTP and OLAP

- Two main categories of database workloads

## Online Transactional Processing (OLTP)

- Point queries and small updates
  <br>

- Read/write a small amount of data (i.e. Record how much of BCS-123 was synthesized, or look up the
  results of an assay ran for BCS-123)
  <br>

  - only need to find BCS-123 and the assay, we don't care about all the other BCS codes
    nor the other assays
    <br>

- high number of short-lived transactions (all the scientists are adding or
  looking up their specific results during the workday)
  <br>

- Postgres, MySQL

## Online Analytical Processing (OLAP)

- Access a large amount or even all of the data (column)
  <br>

- Aggregations, large joins
  <br>
- i.e. What are the most common reactions ran in our labs? What is the average number of steps in the syntheses of our
  molecules?
  <br>

  - need to go through all the reactions or syntheses and count them up
    <br>

- Bulk updates
  <br>

- Data warehouses, data lakes, lakehouses

---

# Very brief intro to OLTP and OLAP

- OLTP and OLAP database systems apply different strategies and techniques to optimize for their specific workload
  <br>

  - B+ Trees, row based storage, tuple-at-a-time query execution
    <br>
  - zone maps, column based storage, vectorized execution
    <br>

- What's good for OLTP, may be a poor choice for OLAP and vice-versa
  <br>
- Wanted to provide context for later on in the talk of why I took a certain approach (Did you try an index?)

---

# duckdb

- in-process (just like importing a library), analytical engine: "the SQLite of analytics"
  <br>

- no client/server:
  <br>

  - client protocol is major bottleneck
    <br>
  - no data transfer or copy costs in many cases
    <br>

- state of the art query execution engine
  <br>

- focused on single-node: "the 99% of problems"
  <br>

- not just limited to it's own data format
  <br>

  - query from MySQL, Postgres, SQLite, remote parquet, csv files, etc...
    <br>

- people have used it many interesting ways besides doing analyses
  <br>

- It's NOT good for every usecase
  <br>

- Gabor Szarnyas - DuckDB (Øredev 2023) https://www.youtube.com/watch?v=6teFN7cwx30
- https://duckdb.org/why_duckdb

---

# duckdb

```python
import duckdb
duckdb.sql("SELECT standard_value, standard_units FROM '../chembl_33_files/activities.parquet' WHERE standard_value < 10 AND standard_units = 'nM'")

┌────────────────┬────────────────┐
│ standard_value │ standard_units │
│     double     │    varchar     │
├────────────────┼────────────────┤
│            6.0 │ nM             │
│            2.0 │ nM             │
│            2.0 │ nM             │
│             ·  │ ·              │
│             ·  │ ·              │
│             ·  │ ·              │
│            6.0 │ nM             │
│            0.5 │ nM             │
│            2.7 │ nM             │
├────────────────┴────────────────┤
│  ? rows (>9999 rows, 6 shown)   │
└─────────────────────────────────┘
```

---

# Cheminformatics workloads and RDKit

- We need software to work with molecules
  <br>
- RDKit: open source cheminformatics software
  <br>
- convert molecular file formats, calculate molecular properties, compare molecules, etc.
  <br>

---

# Implementation of duckdb_rdkit

- Why? Enable cheminformatic workloads in duckdb

---

# Implementation of duckdb_rdkit

##### Registering the `Mol` type in duckdb

```
~~~graph-easy --as=boxart
[ CCO (SMILES) ]  <-> [ RDKit Mol object (in-memory) ]  <-> [ Serialize to binary (on-disk)]
~~~
```

- Need to tell duckdb how to convert this binary into a `Mol` type
- When duckdb sees a `Mol` it can then apply other logic (i.e. can't add two `Mol`)
- Cheminformatics work happens on the in-memory RDKit Mol object

##### Implement basic molecule format functions

- `mol_from_smiles`
- `mol_to_smiles`

---

# Implementing `is_exact_match`

- A basic functionality is to find out if two molecules are the same
- example: SMILES can be written differently, or what if you want to compare a
  molecule represented as SMILES, and another one represented in SDF

```shell
Cc1ccccc1 is toluene

c1ccccc1C is toluene
```

---

# Implementing `is_exact_match`

- Initial implementation: copy the code from the chemicalite RDKit SQLite extension and Postgres RDKit extension
- Example:

  - Query 1:

  ```sql

  select * from molecule where is_exact_match(mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');
  ```

| Query | Standard method (s) |
| :---- | :------------------ |
| 1     | 17.238              |
| 2     | 12.555              |
| 3     | 22.196              |
| 4     | 12.245              |

- default duckdb settings: AMD Ryzen 5 4500U CPU, 16GB RAM, Samsung PM991 SSD
- chembl 33 database ~2.3 million molecules

---

# Implementing `is_exact_match`

##### Initial attempt

- Really poor performance
- What about an index? OLAP uses different indexing strategies than OLTP
  - No B+Tree
- Need to apply different strategies in OLAP

---

# Implementing `is_exact_match`: Umbra Mol

- Umbra-style strings (Neumann, T., Freitag M. CIDR 2020)

```
~~~graph-easy --as=boxart
[length (4B) | prefix (4B) | offset or pointer (8B)]
~~~
```

---

# Implementing `is_exact_match`: Umbra Mol

- Umbra-style strings (Neumann, T., Freitag M. CIDR 2020)

```
~~~graph-easy --as=boxart
[length (4B) | prefix (4B) | offset or pointer (8B)]
~~~
```

- Often times the database system will be comparing strings that don't match
- Don't need to look at the whole string to see 'Hello everyone' and 'Goodbye everyone' are not the same
- Allows for short-circuiting comparisons and speedup
  - pushing 16 bytes through the system instead of thousands of bytes if inlined
    or constantly dereferencing the pointer

---

# Implementing `is_exact_match`: Umbra Mol

- Analyzed the code from Postgres and SQLite RDKit extensions
- The first few checks of the comparison function are short-circuiting checks

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

  // if m1 is substruct of m2 and m2 is substruct of m1, likely to be the same
  // molecule
  RDKit::MatchVectType matchVect;
  bool recursion_possible = false;
  bool do_chiral_match = false;
  bool ss1 = RDKit::SubstructMatch(m1, m2, matchVect, recursion_possible,
                                   do_chiral_match);
  bool ss2 = RDKit::SubstructMatch(m2, m1, matchVect, recursion_possible,
                                   do_chiral_match);
  if (ss1 && !ss2) {
    return false;
  } else if (!ss1 && ss2) {
    return false;
  }

  // the above can still fail in some chirality cases
  std::string smi1 = RDKit::MolToSmiles(m1, do_chiral_match);
  std::string smi2 = RDKit::MolToSmiles(m2, do_chiral_match);
  return smi1 == smi2;
}
```

---

# Implementing `is_exact_match`: Umbra Mol

- Idea: apply the prefix idea from Umbra-style strings
  - pre-compute these "counts" and store it in a prefix
- Umbra-Mol:

```
~~~graph-easy --as=boxart
[Number of atoms | Number of bonds| ... | binary RDKit molecule]
~~~
```

- When comparing two molecules, first compare the prefix, if no match, bail
  out early

```cpp
struct umbra_mol_t {
  uint32_t num_atoms;
  uint32_t num_bonds;
  uint32_t amw;
  uint32_t num_rings;
  uint32_t bmol_size;
  std::string bmol;
}

```

```cpp

bool umbra_mol_cmp(umbra_mol_t m1, umbra_mol_t m2) {
  // check the prefix
  // if any of these values are not equal between the two molecules,
  // there is no way the molecules are the same
  if (m1.num_atoms != m2.num_atoms || m1.num_bonds != m2.num_bonds ||
      m1.amw != m2.amw || m1.num_rings != m2.num_rings) {
    return false;
  }

  // otherwise, run a full check on the molecule objects
  std::unique_ptr<RDKit::ROMol> left_mol(new RDKit::ROMol());
  std::unique_ptr<RDKit::ROMol> right_mol(new RDKit::ROMol());


  // deserialization step
  RDKit::MolPickler::molFromPickle(m1.bmol, *left_mol);
  RDKit::MolPickler::molFromPickle(m2.bmol, *right_mol);
  return mol_cmp(*left_mol, *right_mol);
}

```

---

# Implementing `is_exact_match`: Umbra Mol

##### Umbra Mol results

| Query | Standard method (s) | Umbra-mol (s) | speedup |
| :---- | :------------------ | :------------ | ------: |
| 1     | 17.238              | 0.496         |  34.75x |
| 2     | 12.555              | 0.473         |  26.54x |
| 3     | 22.196              | 0.364         |  60.98x |
| 4     | 12.245              | 0.359         |  34.11x |

- This suggests that the deserialization of the binary molecule to the in-memory
  representation is very expensive

---

# Umbra Mol version 2

- Shared the results on hacker news and got great feedback and ideas from Andrew
  Dalke
- Tried one of his ideas, and found a big regression in performance
- Improved Umbra Mol:
  - shrink the "count" prefix
  - add a filter for bailing out early for substructure matches: dalke_fp
  - store pointer to binary molecule instead of inlining the binary

---

# Umbra Mol version 2

##### improving the "count" prefix

- Umbra Mol version 1 had a 20 byte count prefix -- that's a lot
- Analyzed the chembl 33 dataset and looked at the distribution of the data
- Can capture 99% of the data with far fewer bytes
  - 99% of the molecules have 8 rings or less. That only requires 3 bits, not 32!

```
D select quantile_cont(num_atoms, 0.99), quantile_cont(num_bonds, 0.99), quantile_cont(num_rings, 0.99), quantile_cont(amw::integer, 0.99) from molecule;
┌────────────────────────────────┬────────────────────────────────┬────────────────────────────────┬───────────────────────────────────────────┐
│ quantile_cont(num_atoms, 0.99) │ quantile_cont(num_bonds, 0.99) │ quantile_cont(num_rings, 0.99) │ quantile_cont(CAST(amw AS INTEGER), 0.99) │
│             double             │             double             │             double             │                  double                   │
├────────────────────────────────┼────────────────────────────────┼────────────────────────────────┼───────────────────────────────────────────┤
│                          200.0 │                           37.0 │                            8.0 │                                    1446.0 │
└────────────────────────────────┴────────────────────────────────┴────────────────────────────────┴───────────────────────────────────────────┘
Run Time (s): real 0.190 user 0.182839 sys 0.164021

# and the binary molecule size

chembl_33=#  select percentile_cont(0.99) within group (order by octet_length(mol_to_pkl(rdkit_mol))) from compound_structures;
-[ RECORD 1 ]---+-----
percentile_cont | 1310

Time: 124568.271 ms (02:04.568)
```

---

# Umbra Mol version 2

##### improving the "count" prefix

Using the chembl analysis as a guide, here is how the new counts prefix looks like:

- number of atoms, 7 bits
- number of bonds, 6 bits
- number of rings, 3 bits
- amw, 11 bits

total: 27 bits (~4B as opposed to 20B)

---

# Umbra Mol version 2

##### dalke_fp as a substructure filter

- Substructure matches: is this structure found in that molecule? Does my molecule
  contain a carboxylic acid?

---

# Umbra Mol version 2

##### dalke_fp as a substructure filter

- If N exists in my query molecule (or substructure), but not in the target molecule (the molecule
  I am comparing my query to), then the query molecule cannot be a substructure of the target

- If N exists in the target, but not in the query molecule, it is still possible however
  that there is something else in the query molecule that does match.

  - Query: CCC

  - Target: CCCCCCN

  - CCC is a substructure of CCCCCCN

  - identify molecular features in a dataset that split the data into smaller and smaller bins

```
                        all mols

                has O               no O

         has N     no N        has CC   no CC

```

- The hard part is identifying these features
  <br>
- Once identified: create a fingerprint (1 if present, 0 if not) and use for short-circuiting
  <br>
- 55 bits were found by Greg Landrum to be effective already
  <br>
- I incorporated these 55 bits, `dalke_fp`, into the Umbra Mol prefix
  <br>

- This suggests a possibility of bailing out early. But how to get a meaningful filter?
- Details here:
  - http://www.dalkescientific.com/writings/diary/archive/2012/06/11/optimizing_substructure_keys.html
  - https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02078.html

---

# Umbra Mol version 2

##### dalke_fp as a substructure filter

- Andrew Dalke came up with a clever way to find features in a molecular dataset that
  when splitting the data on these features, the molecules end up in much smaller buckets.

(Sorry: it's a tree but it fell over)

```
~~~graph-easy --as=boxart
[all molecules] -> [has O]
[all molecules] -> [does not have O]
[has O] -> [has N]
[has O] -> [does not have N]
[does not have O] -> [has Cl]
[does not have O] -> [does not have Cl]
[has Cl] -> [...]
[does not have Cl] -> [and so on ...]
~~~
```

- The presence of these features can be marked in a bit vector or fingerprint

- Details here:
  - http://www.dalkescientific.com/writings/diary/archive/2012/06/11/optimizing_substructure_keys.html
  - https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02078.html

---

# Umbra Mol version 2

##### dalke_fp as a substructure filter

- The hardest part is figuring out what are the meaningful features to split on
- Thankfully Andrew Dalke & Greg Landrum of RDKit worked on this idea and shared
  their work
- They found that 55 bits from Andrew Dalke's calculations worked well in bailing
  out early and avoiding a full substructure test
- I copied these features and put it into Umbra Mol

- Details here:
  - http://www.dalkescientific.com/writings/diary/archive/2012/06/11/optimizing_substructure_keys.html
  - https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02078.html

---

# Umbra Mol version 2

##### dalke_fp as a substructure filter

The prefix in UmbraMol now looks like this:

```
~~~graph-easy --as=boxart
[count prefix (4B) | dalke_fp (8B)|... ]
~~~

```

---

# Umbra Mol version 2

##### storing a pointer to binary molecule instead of inlining

- Another key idea of Umbra-style strings is to store a pointer to the string
- Umbra Mol version 1 inlined the binary molecule
- My experiments (not shown) suggested that searches were getting slower because
  the Umbra Mol was getting too big (tried to put in the counts and dalke_fp and
  binary molecule all in the struct)
- Idea: do what Umbra-style strings do instead of inlining the binary molecule

---

# Umbra Mol version 2

##### storing a pointer to binary molecule instead of inlining

- But, this pointer is a virtual memory address
- When the data goes to disk, this memory address is no longer valid.
- Pointer swizzling
  - translate between virtual memory address and, for example, an offset on a
    database page
- How to do pointer swizzling as an extension?
  - I assume it would need to be very integrated into how the data is
    represented on disk and the buffer manager

---

# Umbra Mol version 2

##### storing a pointer to binary molecule instead of inlining

- duckdb implements Umbra-style strings
  - `string_t`
  - that means they must do pointer swizzling
- After a tip from someone in duckdb's discord server, I dug around the code
  and found what kind of data gets pointers swizzling applied to it
- Long story short: if I can make Umbra Mol look like a `string_t`, duckdb will
  handle this heavy lifting for me

---

# Umbra Mol version 2

##### storing a pointer to binary molecule instead of inlining

```
~~~graph-easy --as=boxart
[Umbra Mol ] -> [string_t] -> [duckdb internal pointer swizzling, etc.] -> [disk]
~~~
```

```
~~~graph-easy --as=boxart
[disk] -> [duckdb internal pointer swizzling, etc.] -> [string_t] -> [Umbra Mol] -> [cheminformatics stuff]
~~~
```

---

# Umbra Mol version 2

##### storing a pointer to binary molecule instead of inlining

- Serializing Umbra Mol:
  - just serialize to binary and store it in a `string_t`

```c++

  std::string buffer;
  // total size = prefix size + binary molecule size
  buffer.reserve(total_size);
  // the counts prefix - 4 bytes
  buffer.append(reinterpret_cast<const char *>(&prefix),
                umbra_mol_t::COUNT_PREFIX_BYTES);
  // dalke_fp part of the prefix - 8 btyes
  buffer.append(reinterpret_cast<const char *>(&dalke_fp),
                umbra_mol_t::DALKE_FP_PREFIX_BYTES);
  // the binary molecule from RDKit
  buffer.append(binary_mol);

```

- `string_t` has a constructor that takes `std::string`

---

# Umbra Mol version 2

##### storing a pointer to binary molecule instead of inlining

- Deserializing Umbra Mol:
- In duckdb, the `string_t` stores data in this struct;

```c++
struct {
    uint32_t length;
    char prefix[4];
    char *ptr;
} pointer;
```

- Umbra Mol version 2 looks like this:
- The main difference is the length of the prefix
  - 4 bytes vs 12 bytes (4B counts + 8B dalke_fp)

```c++
  struct {
    uint32_t length;
    char prefix[PREFIX_LENGTH];
    const char *ptr;
  } value;

```

- Idea: just copy the first 12 bytes (PREFIX_LENGTH) of the data in `ptr` in `string_t`
  to Umbra Mol's prefix field, and then set the `ptr` of Umbra Mol to `ptr` in `string_t`

---

# Umbra Mol version 2

##### storing a pointer to binary molecule instead of inlining

- Deserializing Umbra Mol:

```c++
  // This constructor is used to convert a `string_t` to an `umbra_mol_t`:
  umbra_mol_t(string_t buffer) {
    value.length = buffer.GetSize();
    // first 12 bytes are the prefix
    memset(value.prefix, 0, PREFIX_LENGTH);
    memcpy(&value.prefix, buffer.GetData(), PREFIX_LENGTH);
    // string_t::GetData() gets the pointer to the string
    value.ptr = buffer.GetData();

    D_ASSERT(value.ptr == buffer.GetData());
  }

```

---

# Umbra Mol version 2

- Put together:
  - Improved prefix: 4B count prefix + 8B substructure filter -- do more with 12B than
    I did with 20B
  - Store pointer to binary molecule instead of inlining it in the struct
    - don't pass the big binary molecule around if we can answer our questions
      with the prefix instead

---

# Umbra Mol version 2 experiments

##### Exact match

| Query | Standard method (s) | Umbra-mol part 2 (s) | speedup (Umbra-mol vs standard method) | Postgres control (s) |
| :---- | :------------------ | :------------------- | :------------------------------------- | :------------------- |
| 1     | 17.238              | 0.179                | 96x                                    | 0.084                |
| 2     | 12.555              | 0.145                | 87x                                    | 233                  |
| 3     | 13.027              | 0.263                | 50x                                    | 2.47                 |
| 4     | 12.245              | 0.255                | 48x                                    | 6.185                |

##### Substructure match

| Query | Standard method (s) | Umbra-mol part 2 (s) | speedup (Umbra-mol vs standard method) | Postgres control (s) |
| :---- | :------------------ | :------------------- | :------------------------------------- | :------------------- |
| 1     | 23.388              | 0.267                | 88x                                    | 0.741                |
| 2     | 14.094              | 5.93                 | 2x                                     | 98                   |
| 3     | 14.294              | 0.553                | 26x                                    | 12.114               |
| 4     | 13.994              | 6.804                | 2x                                     | 1237 (20 min)        |

- duckdb and postgres were just run on default settings
