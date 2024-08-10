---
layout: post
title: Umbra-style molecules - part 2
description: ""
summary: ""
tags: [chembl, chemistry, rdkit, databases, duckdb, umbra]
---

# Abstract

In part 1, I applied the prefix idea from Umbra-style strings to molecules and improved
exact match queries via short-circuiting in duckdb_rdkit. I then tried to add
more data to the prefix to also enable short-circuiting for substructure matches.
This caused a degradation in exact match performance. Experiments suggested that
execution time increased as the Umbra-mol struct got larger (more data in the prefix
plus the binary molecule). Here in part 2, I improved the prefix by fitting more
useful data into 12 bytes than the 20 bytes I used in the initial implementation.
Furthermore, I applied the second key idea from Umbra-style
strings -- storing a pointer to the full binary molecule, rather than inlining
it in the struct. This allowed for speedup of not only exact matches, but also
enabled substructure matches when combined with a substructure filter developed by dalke.

Here are the results, and I describe the process more below:

Exact match experiments (Standard method, Umbra-mol are in the duckdb_rdkit extension):

| Query | Standard (s) | Umbra-mol (s) | speedup (Umbra-mol vs standard method) | Postgres control (s) |
| :---- | :----------- | :------------ | :------------------------------------- | :------------------- |
| 1     | 17.238       | 0.179         | 96x                                    | 0.084                |
| 2     | 12.555       | 0.145         | 87x                                    | 233                  |
| 3     | 13.027       | 0.263         | 50x                                    | 2.47                 |
| 4     | 12.245       | 0.255         | 48x                                    | 6.185                |

Substructure match experiments (Standard method, Umbra-mol are in the duckdb_rdkit extension):

| Query | Standard (s) | Umbra-mol (s) | speedup (Umbra-mol vs standard method) | Postgres control (s) |
| :---- | :----------- | :------------ | :------------------------------------- | :------------------- |
| 1     | 23.388       | 0.267         | 88x                                    | 0.741                |
| 2     | 14.094       | 5.93          | 2x                                     | 98                   |
| 3     | 14.294       | 0.553         | 26x                                    | 12.114               |
| 4     | 13.994       | 6.804         | 2x                                     | 1237 (20 min)        |

# Introduction

Duckdb is an in-process analytical database system, and I began experimenting with
adding cheminformatic capabilities to duckdb as a duckdb extension. Initial experiments
showed that exact match molecule queries were quite slow.

If you aren't familiar with OLAP and OLTP, I'll try give a brief, brief intro.
Online transactional processing (OLTP) databases are optimized for workloads
that require point lookups. These workloads may just want to access one or
few records in the database, like adding a molecule to the
database, adding some assay result for that molecule, or finding out what assays
were run on that molecule, for example.
Online analytical processing (OLAP) databases are optimized for workloads that
need to touch many or all records, like aggregation queries where you want to
compute an average for a column.

There's a lot to say about the strategies employed by OLAP systems, but I will
just touch on one part because the question may arise in relation to the topic
of this post: "did you try using an index to speed up the searches?"

OLAP database systems like duckdb typically do not use the indexes commonly found in
OLTP databases (e.g. Postgres), like B-Trees. B-Trees are very powerful for getting
in and out fast to find the record you're interested in -- point lookups.
However, B-Trees may not be a good fit for OLAP workloads where the system may need to scan
all the records to find an answer to your query, or for bulk updates where the overhead
of maintaining a B-Tree is too expensive.

There are many interesting strategies OLAP systems employ to make full table
scans really fast. These strategies connect to one another and touch on
many different parts of the database system, from the storage format
(columnar format instead of row-wise) to the query execution engine
(vectorized instead of tuple at a time). I recommend the database systems courses available
on Youtube by the professors at CMU, if you are interested in more information.

That being said, what strategies can be employed to speed up these cheminformatic comparisons,
in duckdb?

I began to look into how to improve this performance, and was inspired by [Umbra-style strings], aka
German-style strings, which is a string format adopted by several data processing
systems, including duckdb. There are two key ideas that I applied from
Umbra-style strings to molecules: a prefix for short-circuting comparisons, and storing
a pointer to a string rather than inlining the string in the in-memory representation.

Molecule comparison is more complicated than strings ([learn more here]), but
the idea of short-circuiting comparisons is not a new one for molecule searches, and
I thought the ideas from Umbra-style strings could be translated to molecules.

In my intial attempt with "Umbra-mols" (see [part 1]), I applied only the prefix idea. I
followed previous examples from [chemicalite] and the [RDKit Postgres extension],
and used counts like the number of atoms, bonds, etc. which are cheap to
compute and store. I put these counts in a prefix, and compared the prefixes of
two molecules to bail out of the function early if the prefixes did not match.
If two molecules have different number of atoms, they cannot be the same molecule.
This gave a huge improvement in exact match comparisons, achieving up to 60x speedup
over the "standard" implementation without a prefix in some simple experiments
I did.

Processing a molecule involves first parsing some kind of string molecular format into
a RDKit molecule object. Then to store this data, it needs to be serialized to binary.
SMILES -> RDKit Molecule object -> serialize to binary (binary molecule)
To do cheminformatic work on this stored molecule, we need to deserialize the
binary into the in-memory RDKit Molecule object.

The prefix sits in front of a binary molecule object, like so:

```
| 4 bytes - number of atoms | 4 bytes - number of bonds | 4 bytes - amw | 4 bytes - number of rings | 4 bytes - size of the binary molecule | n bytes - binary molecule |
```

The initial experiments suggested that this deserialization process is very expensive,
especially in an OLAP scenario where the full column is being scanned. Every
binary molecule would be deserialized in order to do even the inexpensive comparisons.
The reason the prefix and short-circuit mechanism is so powerful is because it allows
the system to save time by avoiding this deserialization. Only when the prefixes match,
deserializing the binary molecule is necessary to do further checks.

After sharing these results, user dalke on Hacker News gave me some great new ideas
to explore. As I worked on these new ideas, it exposed limitations of the first
implementation. More experiments suggested it was the size of Umbra-mol in-memory representation.
And that led me to explore how to implement the second key idea from Umbra-style strings:
store a pointer to the string rather than inlining the string into the struct.

<!-- To summarize my follow up experiments: -->
<!---->
<!-- - First, I was able to squeeze the counts prefix into 4-bytes. -->
<!--   My initial attempt used 20 bytes for the prefix.
<!---->
<!-- - I tried to extend the Umbra-mol prefix to include a small fingerprint that could be used -->
<!--   to short-circuit substructure matches, a substructure filter. I learned that this is -->
<!--   and idea that dalke and others have thought about at least since [2012]. -->
<!--   I used the same [55 bit] fingerprint that dalke calculated ([details here]) -->
<!--   and that Greg at RDKit experimented with, and put that in the Umbra-mol prefix. -->
<!--   However, this started to degrade the exact match performance a lot. -->
<!--   Some queries that executed in < 1 sec were now taking > 5 sec. -->

# Improvements to the first attempt at Umbra-mol

## Prefix improvements

I shrank the counts prefix used for short-circuiting exact match comparisons, and
added `dalke_fp`, a set of fingerprints that dalke calculated, for short-circuiting
substructure matches.

### Counts prefix for short-circuiting exact match comparisons

Initially I just used 4 bytes for each field of the counts properties that are used
for short-circuting the exact match comparisons forming a 20 byte prefix.
I did an analysis of the chembl
database molecules, found that the distribution for some of these properties
were extremely skewed, and that I could capture 99% of the dataset with much smaller
integer representations because the values are not very large [see supplementary - 99th percentile](#99percentile).
For example, 99% of the molecules have less than 8 rings, so I could just
use 3 bits for that.

Using the chembl analysis as a guide, here is how the new counts prefix looks like:

- number of atoms, 7 bits
- number of bonds, 6 bits
- number of rings, 3 bits
- amw, 11 bits

Any molecules that come in with values larger than the allotted
bits can represent are simply capped at the max value for that number of bits.
Any collisions with
other molecules is ok because it's just for short-circuiting comparisons, not creating
cryptographic hashes, and since these bits cover 99% of the values they are
responsible for,
the number of collisions due to being capped are expected to be small.

This results in a counts prefix that fits within 4 bytes.

### dalke_fp as a substructure filter

In order to short-circuit substructure matches, we need a way to rule out quickly
that it is impossible for a query molecule to be a substructure of a target molecule.
dalke has written about this in [2012]. Paraphrasing from dalke's article, this is
the idea:

For example, we can have a bit vector, or molecular fingerprint, 1 for the presence of a chemical motif,
and 0 when absent.

If the fragment exists in the query but not in the target, so 1 in the query, 0
in the target, there is no way for a match.

If the fragment exists in the target, but not in the query, it is still
possible there is something else in the query that matches the target, because the common
substructure is not captured in the fingerprint. An example from dalke's article:
"if "N" does not exist in the query then it might still match a target molecule
which has an "N" in it. For example, the query "O" should still match the target
"N#\[N+\]\[O-\]" (nitrous oxide)."

If all fragments that have on bits in the query are also on in the target,
it does not mean that the query is a substructure. It is possible
that there is something present in the query not captured by the fingerprint,
but not in the target. For example, if the target is "NC", and the query is
"NCCCCCCCC", and the fingerprint only accounts for the presence of "N", both
the query and target would have the bit corresponding to "N" set, but the query
is not a substructure of the target.

Thus, it is only possible to short-circuit in the false case, not in the
true case with this substructure filter.

The hardest part here is calculating a good set of fingerprints that can screen
out a large portion of the dataset. dalke calculated a set of fingerprints in
some experiments years ago, and I just used that. See [dalke's article] for the details of
how this set is calculated.

This set, which I call `dalke_fp` fits in 8 bytes, and I put that in the front
of the binary molecule to create a prefix of 12 bytes -- 4 bytes for the counts,
8 bytes for `dalke_fp`. It might be that the `dalke_fp` is also sufficient for
short-circuiting exact match comparisons, but I haven't checked this. If so,
the counts prefix is unnecessary, and the prefix could be shrunk to 8 bytes.

## Storing a pointer to the binary molecule

The second key idea of Umbra-style strings is to store a pointer to the string
rather than inlining the string. Likewise, in this work I store a pointer to the entire
Umbra Mol which has the binary molecule.

My initial attempt inlined the binary molecule in the struct, and that takes up
a lot of space (most binary molecules are over 1,000 bytes [supplementary figure](#99percentile).)
Then when I tried putting the dalke_fp in the prefix, I started to see degradation
of the exact match performance. Some queries that executed in < 1 sec when
the prefix was only counts, were now taking > 5 sec. I did some experiments
that suggested the data per Umbra Mol was getting too big. Because the binary
molecule is inlined, it is being pushed around the system even though it is
not always needed if the prefix enables a short-circuit.

In Umbra-style strings they store a pointer to the string, if it is a "long" string.
Here is the struct of the `string_t` type in duckdb's string implementation.

```c++

	union {
		struct {
			uint32_t length;
			char prefix[4];
			char *ptr;
		} pointer;
		struct {
			uint32_t length;
			char inlined[12];
		} inlined;
	} value;

```

If the pointer is stored, then there is a further step of pointer swizzling that needs
to take place. We need a way to translate that pointer
to a memory address to something that can point to an offset somewhere on disk when
the data moves from memory to disk, and vice versa. I assume that this process is
probably very closely tied to
rest of the duckdb internals like the buffer pool and storage format, and I wasn't sure
if I could alter those internals from an extension in order to do pointer swizzling.

I wanted to see how duckdb did this for the `string_t` type, and I asked the
duckdb Discord channel to see if someone could give me a suggestion
of where to look. I quickly got some help, and that led me to `column_data_collection_segment.cpp`
and a `UnswizzlePointers` function.

I found that this function gets called if the internal type of the data is
a `PhysicalType::VARCHAR`. In the extension, I create `LogicalTypes`, and so
I checked to see which logical types map to a `PhysicalType::VARCHAR`. In `types.cpp`,
I found that `LogicalTypeId::BLOB` maps to `PhysicalType::VARCHAR`, which is already
the logical type I was using for the Umbra-mol.

With this in mind, my idea was that if I can translate my Umbra-mol to a `string_t`,
then the heavy lifting would be done by duckdb. After all it's all just binary
data. The main difference between the Umbra-mol, or `umbra_mol_t` as the struct
in my code is called, and `string_t` is the size of the prefix; I want to store
more data in the prefix than the 4 bytes in `string_t`

Converting an Umbra-mol to `string_t` wasn't really much of a conversion at all.
I simply made the calculations for the prefixes, and packed it into a `std::string` along with the
binary molecule at the end, and then passed that to a `string_t`
constructor.

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

Duckdb's `string_t` implementation does the rest. When it goes
to and from disk, all the swizzling, unswizzling, and other complicated tasks are
taken care of by duckdb since it just sees a `string_t`.

Converting from a `string_t` to `umbra_mol_t` was also rather straightforward,
after I figured out all my pointer reference errors:

Here's the `value` struct inside the `umbra_mol_t` struct:

```c++
  struct {
    uint32_t length;
    char prefix[PREFIX_LENGTH];
    const char *ptr;
  } value;

```

This constructor is used to convert a `string_t` to an `umbra_mol_t`:

```c++

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

Just copy the length (this is the prefix + binary molecule), and then copy
the first 12 bytes (more on the prefix in the next section) into the `prefix`
field, and then simply set the pointer in `umbra_mol_t` to the pointer in
`string_t` that points to the binary molecule.

To get the binary molecule for when further calculations need to be made because
the prefix cannot short-circuit the comparisons, I can copy the bytes that
follow the prefix.

```c++

  std::string GetBinaryMol() {
    idx_t bmol_size = value.length - PREFIX_LENGTH;
    std::string buffer;
    buffer.resize(bmol_size);
    if (value.ptr && value.length > PREFIX_LENGTH) {
      memcpy(&buffer[0], &value.ptr[PREFIX_LENGTH], bmol_size);
    }
    return buffer;
  }

```

The process is then:
binary data string with the data of
an Umbra-mol -> `string_t` -> duckdb internals -> `string_t` -> `umbra_mol_t`

Now with all the pieces in place, let's see how exact and substructure matches
perform now that the binary molecule is not inlined into the `umbra_mol_t`.

# Experiments:

I made up some queries, and for substructure queries, I chose structures from a
list of substructure queries gathered from real queries that dalke has
shared [here](https://hg.sr.ht/~dalke/sqc/browse/README?rev=tip).

## Exact match

Queries are shown in [supplementary 1](#supplementary1).

| Query | Standard method (s) | Umbra-mol part 2 (s) | speedup (Umbra-mol vs standard method) | Postgres control (s) |
| :---- | :------------------ | :------------------- | :------------------------------------- | :------------------- |
| 1     | 17.238              | 0.179                | 96x                                    | 0.084                |
| 2     | 12.555              | 0.145                | 87x                                    | 233                  |
| 3     | 13.027              | 0.263                | 50x                                    | 2.47                 |
| 4     | 12.245              | 0.255                | 48x                                    | 6.185                |

## Substructure match

Queries are shown in [supplementary 2](#supplementary2).

| Query | Standard method (s) | Umbra-mol part 2 (s) | speedup (Umbra-mol vs standard method) | Postgres control (s) |
| :---- | :------------------ | :------------------- | :------------------------------------- | :------------------- |
| 1     | 23.388              | 0.267                | 88x                                    | 0.741                |
| 2     | 14.094              | 5.93                 | 2x                                     | 98                   |
| 3     | 14.294              | 0.553                | 26x                                    | 12.114               |
| 4     | 13.994              | 6.804                | 2x                                     | 1237 (20 min)        |

Note on query 4: I am detecting more matches (16,543 vs 16,165) in the duckdb queries
than with the Postgres RDKit extension.
I get the same values for the standard method and the Umbra-mol which suggests
that the prefix is not causing the behavior. It might be that I have some different
parameters for the `RDKit::SubstructMatch` function, but I tried with the same
parameters as I found in the Postgres RDKit extension code, as well as with chemicalite
(I am following the chemicalite parameters for now), and that didn't solve the
discrepency. I don't know if these are false positives, or if the Postgres RDKit extension
gets false negatives. I have no idea what is the reason,
and I only observe this behavior for this query so far. All the other queries I tried
did not exhibit this behavior.

# Conclusion

Storing the pointer to the entire binary molecule rather than inlining the
data into the `umbra_mol_t` struct has a massive effect on performance.
It also enables the Umbra-mol to pack more useful data into the prefix to
also improve substructure searches.

My initial attempt, kept an inlined binary molecule in the struct.
I was seeing > 10 sec execution times on Substructure match Queries 1 and 2 even
with the `dalke_fp` in the prefix for short-circuiting.
Storing a pointer to the binary molecule drops those down to 0.267 and ~6 seconds
respectively! For exact matches, I saw improvements of ~26-60x on exact match queries over the
standard, non-prefixed, implementation. With the pointer to the binary molecule,
I see 48x-96x speedups, all while having more useful information in the prefix.

The `dalke_fp` seems very effective for short-circuiting, but it also depends on
the molecules. In the limited samples I tried, I found for some substructures the
`dalke_fp` could short-circuit all but ~400 molecules out of ~2.7 million molecules.
That means the vast majority of that search did not even need to dereference the pointer
to the binary molecule nor deserialize it.

In the substructure matches, the query structure has a big effect on the execution time;
the more general the query structure is, the more that it will match molecules in
the database. Since there are many matches, the full deserialization and substructure
match needs to be done to confirm the match. The rarer the query substructure is, the more
the prefix will enable the function to bail out early, and the faster the query will
be.

# Supplementary figures

## <a name="99percentile"></a> 99th percetile of the features used in the count prefix

It's possible to represent most of this data with just a few bits in some cases.

```sql
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

## <a name="supplementary1"></a> Exact match queries

Query 1

```sql

-- umbra mol
select molregno,canonical_smiles, rdkit_mol, umbra_mol from molecule where umbra_is_exact_match(umbra_mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');

-- standard
select molregno,canonical_smiles, rdkit_mol, umbra_mol from molecule where is_exact_match(rdkit_mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');

-- postgres
select molregno,canonical_smiles, rdkit_mol from compound_structures where rdkit_mol@='Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O';

```

Query 2

```sql
-- umbra mol
select molregno,canonical_smiles, rdkit_mol, umbra_mol from molecule where umbra_is_exact_match(umbra_mol,'CCC');
-- standard
select molregno,canonical_smiles, rdkit_mol, umbra_mol from molecule where is_exact_match(rdkit_mol,'CCC');

-- postgres
select molregno,canonical_smiles, rdkit_mol from compound_structures where rdkit_mol@='CCC';

```

Query 3

```sql

-- umbra mol
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.umbra_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE umbra_is_exact_match(m.umbra_mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');

-- standard
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.rdkit_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE is_exact_match(m.rdkit_mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');

-- postgres
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.rdkit_mol FROM compound_structures m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE m.rdkit_mol@='COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC';

```

Query 4

```sql

-- umbra mol
  SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.umbra_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE umbra_is_exact_match(m.umbra_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
      GROUP BY m.umbra_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

-- standard
SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE is_exact_match(m.rdkit_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;


-- postgres
SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM compound_structures m
  INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
  WHERE m.rdkit_mol@='CC(=O)Nc1nnc(S(N)(=O)=O)s1'
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

```

## <a name="supplementary2"></a> Substructure match queries

Query 1

```sql

-- umbra mol
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.umbra_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE umbra_is_substruct(m.umbra_mol, 'COc1cc2c(Nc3cc(CC(=O)Nc4cccc(F)c4F)[nH]n3)ncnc2cc1OCCCN(CCO)CC(C)C');

-- standard
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.rdkit_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE is_substruct(m.rdkit_mol, 'COc1cc2c(Nc3cc(CC(=O)Nc4cccc(F)c4F)[nH]n3)ncnc2cc1OCCCN(CCO)CC(C)C');

-- postgres
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.rdkit_mol FROM compound_structures m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE m.rdkit_mol@='COc1cc2c(Nc3cc(CC(=O)Nc4cccc(F)c4F)[nH]n3)ncnc2cc1OCCCN(CCO)CC(C)C';

```

Query 2

```sql
-- umbra mol
SELECT count(*) FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      WHERE umbra_is_substruct(m.umbra_mol, 'O=CNCCc1ccccc1');

-- standard
SELECT count(*) FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      WHERE is_substruct(m.rdkit_mol, 'O=CNCCc1ccccc1');

-- postgres
SELECT count(*) FROM compound_structures m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      WHERE m.rdkit_mol@>'O=CNCCc1ccccc1';

```

Query 3

```sql

-- umbra mol
SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.umbra_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE umbra_is_substruct(m.umbra_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
      GROUP BY m.umbra_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

-- normal
SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE is_substruct(m.rdkit_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;


-- postgres
SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM compound_structures m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE m.rdkit_mol@>'CC(=O)Nc1nnc(S(N)(=O)=O)s1'
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

```

Query 4

```sql

-- umbra mol
  SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.umbra_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE umbra_is_substruct(m.umbra_mol, 'N1C=CC=N1')
      GROUP BY m.umbra_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

-- standard
  SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE is_substruct(m.rdkit_mol, 'N1C=CC=N1')
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;


-- postgres
SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM compound_structures m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE m.rdkit_mol@>'N1C=CC=N1'
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

```

[Umbra-style strings]: https://cedardb.com/blog/german_strings/
[learn more here]: https://depth-first.com/articles/2019/04/12/the-smiles-substructure-search-fallacy/
[part 1]: https://bodowd.github.io/2024/07/23/umbra-style-molecules
[chemicalite]: https://github.com/rvianello/chemicalite
[RDKit]: https://www.rdkit.org/docs/Overview.html
[RDKit Postgres extension]: https://www.rdkit.org/docs/Cartridge.html
[2012]: http://www.dalkescientific.com/writings/diary/archive/2012/06/11/optimizing_substructure_keys.html
[55 bit]: https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02078.html
[details here]: http://www.dalkescientific.com/writings/diary/archive/2012/06/11/optimizing_substructure_keys.html
[dalke's article]: http://www.dalkescientific.com/writings/diary/archive/2012/06/11/optimizing_substructure_keys.html
