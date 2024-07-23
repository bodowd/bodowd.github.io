---
layout: post
title: making exact molecule match faster in duckdb inspired by umbra
description: ""
summary: ""
tags: [chembl, chemistry, rdkit, databases, duckdb, umbra]
---

# Moleküle? Umbra-style strings applied to molecules to speed up molecular equality searches in duckdb_rdkit

The current implementation of `is_exact_match` for finding molecules in duckdb_rdkit
uses the standard molecule comparision algorithm found in the RDKit Postgres
extension and the RDkit SQLite extension, chemicalite[1]. I was inspired
by Umbra-style strings (or German-style strings[2]) and applied that idea to
the exact match for molecules function in order to speed up the search.

### Umbra-style/German-style strings

Umbra-style or German-style strings[2] is a string storage format that was introduced by Umbra[3]
which makes searching through strings in database systems faster, and this format
has been adopted in several data processing systems.

The Umbra-style string has one representation for short strings,
and one representation for long strings. You can see the details in [2] and [3].
One strategy the Umbra-style string uses is to put a short prefix of the string
in the encoding, and then store a pointer to the rest of the string (in the long string case).

By making the prefix easily accessed, the database system can rapidly rule out
records that have no chance of being equal. If one string starts with "Hi"
and the other string starts with "Dude!", it doesn't matter what the rest of
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
of inexpensive checks to start ruling out if two molecules could be an exact match -- a
fail fast approach.
For example, do the two molecules have the same number of atoms? If not, there
is no way the two molecules are the same. Only after the inexpensive checks are
done, more expensive checks like substructure searches and converting the molecule
to a RDKit canonicalized SMILES string are done.

The molecules are stored as a binary format in the database, and then they need to
be deserialized into RDKit molecule objects in order to access data like the
number of atoms in the molecule, through the molecule object.

This first part of the equality check sounds exactly like the approach taken in
Umbra-style strings, and so I wondered how the exact match would perform if I
pre-calculate these simple values when the structures are loaded into the database,
and then store these values alongside the serialized RDKit molecule object.
This would be similar to the prefix idea of the Umbra-style strings.

Although I did quick checks of how much time the deserialization process takes,
I don't have any data to present which would compare the time for deserialization relative
to the other operations carried out by duckdb in query execution. That would be
interesting to quantify though.

Anyhow, my hypothesis was that by precalculating and storing these values, I could speed
up the search, assuming that the deserialization process is "slow" relative to the
other operations that are carried out during the query execution, by just comparing
these values in a prefix without deserializing the molecule object -- basically
applying the same idea from Umbra-style strings.

### Approach

To implement this idea, I created a new struct to store the prefix, or header,
containing the values that will be used for the "fail-fast" checks in front
of the binary molecule object. Excerpts of the code are shown instead of all
the code to be brief. If you are interested in seeing more, please see the respository[4].

Here is the key idea -- the first 20 bytes of the Umbra-style molecule stores
the pre-calculated values

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

Then, to check if two molecules are an exact match, I can check the pre-calculated
values.

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

That's the whole idea right there.

The additional code to implement this as a duckdb extension can be found in the
repository[4], which I leave out to be brief.

### Experiments

To run the experiment, I first copied the `molecule` table from the chembl_33
Postgres database to a duckdb file using the Postgres duckdb extension.
That contains ~2.4 million records.

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

```cpp

```

#### Query 1:

Using the standard `is_exact_match` algorithm, it takes ~17 seconds. The search
is insensitive to chirality, so multiple hits are returned.

```sql

select * from molecule where is_exact_match(mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');

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

#### Query 2:

Using the standard algorithm, it takes ~12 seconds.

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

[1]: https://github.com/rvianello/chemicalite
[2]: https://cedardb.com/blog/german_strings/
[3]: https://db.in.tum.de/~freitag/papers/p29-neumann-cidr20.pdf
[4]: https://github.com/bodowd/duckdb_rdkit/tree/umbra-style-mol
