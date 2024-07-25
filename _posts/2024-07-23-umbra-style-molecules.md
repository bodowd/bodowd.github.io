---
layout: post
title: Umbra-style/German-style strings applied to molecules speed up exact match queries on molecules in duckdb_rdkit
description: ""
summary: ""
tags: [chembl, chemistry, rdkit, databases, duckdb, umbra]
---

The current implementation of `is_exact_match` for finding molecules in duckdb_rdkit
uses the standard molecule comparision algorithm found in the RDKit Postgres
extension and the RDkit SQLite extension, chemicalite [[1]].
Here, I apply ideas from Umbra-style/German-style strings [[2]] to speed up
exact search on molecules in duckdb by ~26-60x depending on the type of query.

Jump to results: [results](#results)

Code: [code]

Hacker news discussion: [hn]

- Thanks to dalke, who has given me helpful feedback and new ideas to look into

### Umbra-style/German-style strings

Umbra-style, or German-style strings, [[2]] is a string storage format that was introduced by Umbra [[3]]
which makes searching through strings in database systems faster, and this format
has been adopted in several data processing systems.

The Umbra-style string has one representation for short strings,
and one representation for long strings. You can see the details in [[2]] and [[3]].
One strategy the Umbra-style string uses is to put a short prefix of the string
in the encoding, and then store a pointer to the rest of the string (in the long string case).

The prefix allows the database system to rapidly rule out
records that have no chance of being equal without following the pointer to the
entire string. If one string starts with "Hi"
and the other string starts with "Bye", it doesn't matter what the rest of
the string is, these two strings are not equal.

I wanted to apply this prefix idea to the exact molecule search.

### Analysis of the exact match algorithm drew a connection to Umbra-style strings

Let's start by looking at the exact match algorithm from RDKit Postgres and chemicalite.
Below, you can see the code that I adapted from chemicalite.

```cpp
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

That algorithm takes two RDKit molecule objects, and then runs a series
of inexpensive checks to start ruling out molecules that cannot be a match and short-circuit
before running more expensive checks.
For example, do the two molecules have the same number of atoms? If not, there
is no way the two molecules are the same. Only after the inexpensive checks are
done, more expensive checks like substructure searches and converting the molecule
to a RDKit canonicalized SMILES string are done.

The molecules are stored in a binary format in the database, and they need to
be deserialized into RDKit molecule objects in order to access data like the
number of atoms in the molecule.

I wondered how the exact match would perform if I
pre-calculate and store these simple values when the structures are loaded into the database,
which would allow short-circuiting without the deserialization of the binary molecule.
This would be similar to the prefix idea of the Umbra-style strings.
Assuming that the deserialization process is "slow" relative to the
other operations that are carried out during the query execution, I thought this could
provide a speedup.

Although I did quick checks of how much time the deserialization process takes,
I don't have any data to present which would compare the time deserialization takes relative
to the other operations carried out by duckdb during query execution. That would be
interesting to quantify though.

### Approach

To implement this idea, I created a new struct to store the prefix, or header,
containing the values that will be used for the quick checks in front
of the binary molecule object. Excerpts of the code are shown here for brevity.
If you are interested in seeing more, please see the respository [[4]].

Here is the key idea -- once serialized the first 20 bytes of the Umbra-style molecule stores
the pre-calculated values, and the rest of the binary molecule is stored afterwards,
in order to be available for the more expensive checks, if needed.

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

When serialized to binary, it becomes something like:

`02 00 00 00 01 00 00 00 1e 00 00 00 00 00 00 00 41 00 00 00 ef be ad de 00 00 00 00 0f 00 00 00 00 00 00 00 00 00 00 00 02 00 00 00 01 00 00 00 80 01 06 00 60 00 00 00 01 03 06 00 60 00 00 00 01 03 0b 00 01 00 14 00 00 00 00 17 04 00 00 00 00 00 00 00 16`

for 'CC'. You can see the first four bytes are `02 00 00 00` representing two atoms,
the next four bytes are `01 00 00 00` representing one bond in the molecule,
etc. Then the sequence starting with `ef be ad de ...` is the binary format of
the RDKit molecule.

This is what will be stored in the database in a `UmbraMol` column in the duckdb_rdkit
extension.

Then, to check if two molecules are an exact match, simply deserialize
the `umbra_mol_t` and check the pre-calculated values.
Only when the fast checks do not confirm that the molecules are unequal (meaning
it is still possible the two molecules are the same),
we can deserialize the binary molecule for the more expensive checks.

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

  RDKit::MolPickler::molFromPickle(m1.bmol, *left_mol);
  RDKit::MolPickler::molFromPickle(m2.bmol, *right_mol);
  return mol_cmp(*left_mol, *right_mol);
}

```

In this way, the deserialization of the binary
RDKit molecule can be done only if it is necessary.
If most of the time there are not duplicate molecules,
it would be unnecessary to deserialize the molecule
object, and therefore hopefully this strategy reduces the execution time of the query.

The additional code (like creating a new type, casts, and functions) required
to implement this as a duckdb extension can be found in the
repository [[4]].

### Experiments

Next, I wanted to see how it performs.

To run the experiment, I first copied the `molecule` table from the chembl_33
Postgres database to a duckdb file using the Postgres duckdb extension.
This table contains ~2.4 million records.

As a control experiment, I also ran the same query in Postgres with the RDKit
extension (after constructing an index on
the RDKit molecule object column).
Note: the column names and table names may be different between the Postgres
examples and the duckdb examples because I
renamed the `compound_structures` table to `molecule` when importing the data
to duckdb.

```sql
> select count(*) from molecule;

┌──────────────┐
│ count_star() │
│    int64     │
├──────────────┤
│      2372674 │
└──────────────┘
```

Next, I created a column type for the normal RDKit molecule type, as well as for
the Umbra-style molecule type that I implemented and populated those columns.

```sql
ALTER TABLE molecule ADD COLUMN mol Mol;
UPDATE molecule SET mol=mol_from_smiles(smiles);

ALTER TABLE molecule ADD COLUMN umbra_mol UmbraMol;
UPDATE molecule SET umbra_mol=umbra_mol_from_smiles(mol);

```

#### Query 1:

Using the standard storage format and `is_exact_match` algorithm, it takes ~17 seconds. The search
is insensitive to chirality, so multiple hits are returned.

```sql

D select * from molecule where is_exact_match(mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');

100% ▕████████████████████████████████████████████████████████████▏
┌─────────┬──────────────────────┬──────────────────────┬─────────────────────────────────────┐
│   id    │        smiles        │         mol          │              umbra_mol              │
│  int32  │       varchar        │         mol          │              umbramol               │
├─────────┼──────────────────────┼──────────────────────┼─────────────────────────────────────┤
│   37202 │ Cc1cn([C@H]2C[C@@H…  │ Cc1cn([C@H]2C[C@@H…  │ Cc1cn([C@H]2C[C@@H](N=[N+]=[N-])[…  │
│  298564 │ Cc1cn(C2CC(N=[N+]=…  │ Cc1cn(C2CC(N=[N+]=…  │ Cc1cn(C2CC(N=[N+]=[N-])C(CO)O2)c(…  │
│  372432 │ Cc1cn([C@H]2C[C@H]…  │ Cc1cn([C@H]2C[C@H]…  │ Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C…  │
│  703067 │ Cc1cn([C@@H]2C[C@H…  │ Cc1cn([C@@H]2C[C@H…  │ Cc1cn([C@@H]2C[C@H](N=[N+]=[N-])[…  │
│ 1825987 │ Cc1cn([C@@H]2C[C@H…  │ Cc1cn([C@@H]2C[C@H…  │ Cc1cn([C@@H]2C[C@H](N=[N+]=[N-])[…  │
│ 2542096 │ Cc1cn(C2C[C@H](N=[…  │ Cc1cn(C2C[C@H](N=[…  │ Cc1cn(C2C[C@H](N=[N+]=[N-])[C@@H]…  │
│   27307 │ Cc1cn([C@H]2C[C@H]…  │ Cc1cn([C@H]2C[C@H]…  │ Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C…  │
└─────────┴──────────────────────┴──────────────────────┴─────────────────────────────────────┘
Run Time (s): real 17.238 user 87.310449 sys 0.688456
```

Using the Umbra-style molecule, it takes 0.496 seconds!

```sql
D select * from molecule where umbra_is_exact_match(umbra_mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');
┌─────────┬──────────────────────┬──────────────────────┬─────────────────────────────────────┐
│   id    │        smiles        │         mol          │              umbra_mol              │
│  int32  │       varchar        │         mol          │              umbramol               │
├─────────┼──────────────────────┼──────────────────────┼─────────────────────────────────────┤
│   37202 │ Cc1cn([C@H]2C[C@@H…  │ Cc1cn([C@H]2C[C@@H…  │ Cc1cn([C@H]2C[C@@H](N=[N+]=[N-])[…  │
│  298564 │ Cc1cn(C2CC(N=[N+]=…  │ Cc1cn(C2CC(N=[N+]=…  │ Cc1cn(C2CC(N=[N+]=[N-])C(CO)O2)c(…  │
│  372431 │ Cc1cn([C@H]2C[C@H]…  │ Cc1cn([C@H]2C[C@H]…  │ Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C…  │
│  703067 │ Cc1cn([C@@H]2C[C@H…  │ Cc1cn([C@@H]2C[C@H…  │ Cc1cn([C@@H]2C[C@H](N=[N+]=[N-])[…  │
│ 1825987 │ Cc1cn([C@@H]2C[C@H…  │ Cc1cn([C@@H]2C[C@H…  │ Cc1cn([C@@H]2C[C@H](N=[N+]=[N-])[…  │
│ 2542096 │ Cc1cn(C2C[C@H](N=[…  │ Cc1cn(C2C[C@H](N=[…  │ Cc1cn(C2C[C@H](N=[N+]=[N-])[C@@H]…  │
│   27307 │ Cc1cn([C@H]2C[C@H]…  │ Cc1cn([C@H]2C[C@H]…  │ Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C…  │
└─────────┴──────────────────────┴──────────────────────┴─────────────────────────────────────┘
Run Time (s): real 0.496 user 2.407703 sys 0.028241
```

Postgres control (with gist index on molecule column)

```sql
chembl_33=# select molregno, canonical_smiles from compound_structures where rdkit_mol@='Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O';
 molregno |                      canonical_smiles
----------+------------------------------------------------------------
    37202 | Cc1cn([C@H]2C[C@@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O
    27307 | Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O
   298564 | Cc1cn(C2CC(N=[N+]=[N-])C(CO)O2)c(=O)[nH]c1=O
  1825987 | Cc1cn([C@@H]2C[C@H](N=[N+]=[N-])[C@H](CO)O2)c(=O)[nH]c1=O
   372431 | Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@H](CO)O2)c(=O)[nH]c1=O
   703067 | Cc1cn([C@@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O
  2542096 | Cc1cn(C2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O
(7 rows)

Time: 84.304 ms
```

#### Query 2:

Using the standard storage format and algorithm, it takes ~12 seconds.

```sql
D select * from molecule where is_exact_match(mol,'CCC');
100% ▕████████████████████████████████████████████████████████████▏
┌────────┬─────────┬─────┬───────────┐
│   id   │ smiles  │ mol │ umbra_mol │
│ int32  │ varchar │ mol │ umbramol  │
├────────┼─────────┼─────┼───────────┤
│ 222072 │ CCC     │ CCC │ CCC       │
└────────┴─────────┴─────┴───────────┘
Run Time (s): real 12.555 user 67.567103 sys 0.090192
```

And 0.473 seconds with the Umbra-style molecule.

```sql
D select * from molecule where umbra_is_exact_match(umbra_mol,'CCC');
┌────────┬─────────┬─────┬───────────┐
│   id   │ smiles  │ mol │ umbra_mol │
│ int32  │ varchar │ mol │ umbramol  │
├────────┼─────────┼─────┼───────────┤
│ 222072 │ CCC     │ CCC │ CCC       │
└────────┴─────────┴─────┴───────────┘
Run Time (s): real 0.473 user 2.480841 sys 0.099796
```

The Postgres control query is really slow...
I'm not sure why

```sql
chembl_33=# explain analyze select molregno, canonical_smiles from compound_structures where rdkit_mol@='CCC';
                                                              QUERY PLAN
--------------------------------------------------------------------------------------------------------------------------------------
 Index Scan using molidx on compound_structures  (cost=0.41..8.43 rows=1 width=63) (actual time=70015.058..232891.019 rows=1 loops=1)
   Index Cond: (rdkit_mol @= 'CCC'::mol)
   Rows Removed by Index Recheck: 1846471
 Planning Time: 0.719 ms
 Execution Time: 232893.433 ms
(5 rows)

Time: 232944.928 ms (03:52.945)
```

## Queries with joining and aggregation

For these queries, I copied a few more tables from the chembl_33 postgres db into duckdb
so that I can test how the performance is when joining tables.

### Query 3:

With the standard storage format and algorithm, it takes ~22 seconds.

```sql
D SELECT pbd.prediction_method, a.value, a.relation, m.mol FROM molecule m
    INNER JOIN activities a ON a.molregno=m.id
    INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
    INNER JOIN compound_properties cp ON cp.molregno=m.id
    WHERE is_exact_match(m.mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');
100% ▕████████████████████████████████████████████████████████████▏
┌───────────────────┬────────┬──────────┬───────────────────────────────────────────────────────────────────────┐
│ prediction_method │ value  │ relation │                                  mol                                  │
│      varchar      │ double │ varchar  │                                  mol                                  │
├───────────────────┼────────┼──────────┼───────────────────────────────────────────────────────────────────────┤
│ Multi domain      │   0.52 │ =        │ COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC │
└───────────────────┴────────┴──────────┴───────────────────────────────────────────────────────────────────────┘
Run Time (s): real 22.196 user 114.707474 sys 1.247569
```

With the Umbra-mol, it takes 0.364 seconds.

```sql
D SELECT pbd.prediction_method, a.value, a.relation, m.umbra_mol FROM molecule m
    INNER JOIN activities a ON a.molregno=m.id
    INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
    INNER JOIN compound_properties cp ON cp.molregno=m.id
    WHERE umbra_is_exact_match(m.umbra_mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');
┌───────────────────┬────────┬──────────┬───────────────────────────────────────────────────────────────────────┐
│ prediction_method │ value  │ relation │                               umbra_mol                               │
│      varchar      │ double │ varchar  │                               umbramol                                │
├───────────────────┼────────┼──────────┼───────────────────────────────────────────────────────────────────────┤
│ Multi domain      │   0.52 │ =        │ COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC │
└───────────────────┴────────┴──────────┴───────────────────────────────────────────────────────────────────────┘
Run Time (s): real 0.364 user 1.994158 sys 0.060483
```

Postgres control:

```sql
chembl_33=# SELECT pbd.prediction_method, a.value, a.relation, m.rdkit_mol FROM compound_structures m
    INNER JOIN activities a ON a.molregno=m.molregno
    INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
    INNER JOIN compound_properties cp ON cp.molregno=m.molregno
    WHERE m.rdkit_mol@='COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC';
 prediction_method | value | relation |                               rdkit_mol
-------------------+-------+----------+-----------------------------------------------------------------------
 Multi domain      |  0.52 | =        | COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC
(1 row)

Time: 161.883 ms
```

### Query 4:

Using the standard storage format and algorithm, it takes ~12 seconds

```sql
D SELECT avg(a.value), count(a.value), a.relation, m.mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.id
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  INNER JOIN compound_properties cp ON cp.molregno=m.id
  WHERE is_exact_match(m.mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
  GROUP BY m.mol, a.relation;
100% ▕████████████████████████████████████████████████████████████▏
┌────────────────────┬──────────────────┬──────────┬────────────────────────────┐
│   avg(a."value")   │ count(a."value") │ relation │            mol             │
│       double       │      int64       │ varchar  │            mol             │
├────────────────────┼──────────────────┼──────────┼────────────────────────────┤
│ 1630.1423282178293 │             1818 │ =        │ CC(=O)Nc1nnc(S(N)(=O)=O)s1 │
└────────────────────┴──────────────────┴──────────┴────────────────────────────┘
Run Time (s): real 12.245 user 64.312706 sys 1.056914
```

Using the Umbra-style storage format and algorithm, it takes ~0.36 seconds.

```sql
D SELECT avg(a.value), count(a.value), a.relation, m.umbra_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.id
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  INNER JOIN compound_properties cp ON cp.molregno=m.id
  WHERE umbra_is_exact_match(m.umbra_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
  GROUP BY m.umbra_mol, a.relation;
┌───────────────────┬──────────────────┬──────────┬────────────────────────────┐
│  avg(a."value")   │ count(a."value") │ relation │         umbra_mol          │
│      double       │      int64       │ varchar  │          umbramol          │
├───────────────────┼──────────────────┼──────────┼────────────────────────────┤
│ 1630.142328217826 │             1818 │ =        │ CC(=O)Nc1nnc(S(N)(=O)=O)s1 │
└───────────────────┴──────────────────┴──────────┴────────────────────────────┘
Run Time (s): real 0.359 user 1.912700 sys 0.077231

```

Postgres control:

```sql

chembl_33=# SELECT avg(a.value), count(a.value), a.relation, m.rdkit_mol FROM compound_structures m
INNER JOIN activities a ON a.molregno=m.molregno
INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
INNER JOIN compound_properties cp ON cp.molregno=m.molregno
WHERE m.rdkit_mol@='CC(=O)Nc1nnc(S(N)(=O)=O)s1'
GROUP BY m.rdkit_mol, a.relation;
          avg           | count | relation |         rdkit_mol
------------------------+-------+----------+----------------------------
 1630.14232821782178218 |  1818 | =        | CC(=O)Nc1nnc(S(N)(=O)=O)s1
(1 row)

Time: 900.934 ms
```

## <a name="results"></a>Results

The results are displayed together in the table below. The `real` run time
is displayed in seconds for the standard method and the Umbra-mol method.
Speedup is calculated by `standard method (s) / Umbra-mol (s)`

| Query | Standard method (s) | Umbra-mol (s) | speedup |
| :---- | :------------------ | :------------ | ------: |
| 1     | 17.238              | 0.496         |  34.75x |
| 2     | 12.555              | 0.473         |  26.54x |
| 3     | 22.196              | 0.364         |  60.98x |
| 4     | 12.245              | 0.359         |  34.11x |

[1]: https://github.com/rvianello/chemicalite
[2]: https://cedardb.com/blog/german_strings/
[3]: https://db.in.tum.de/~freitag/papers/p29-neumann-cidr20.pdf
[4]: https://github.com/bodowd/duckdb_rdkit/tree/umbra-style-mol
[results]: https://github.com/bodowd/umbra-style-molecules
[code]: https://github.com/bodowd/duckdb_rdkit/tree/umbra-style-mol
[hn]: https://news.ycombinator.com/item?id=41061014
