---
layout: post
title: rdkit extension for duckdb
description: ""
summary: ""
tags: [duckdb, rdkit, databases]
---

I have been trying to build an extension ([duckdb_rdkit]) for [duckdb] that gives it the ability
to do cheminformatic work via [RDKit].

I've been following and learning from the duckdb extension [spatial],
the sqlite3 extension [chemicalite], as well as the RDKit Postgres extension
as a pattern for building the extension. Many thanks to their work.

## Working with molecules in a computer

Working with molecules requires some special handling. Here are some nice resources
that go into more details: [Daylight], or [depth first blog].

As a quick example, a molecule is a three dimensional dynamic structure in reality
but we often draw it on paper like a 2D graph with nodes being the atoms and
edges being the bonds.

One way to represent a molecule is a SMILES string, a string representation of
a molecule built from a depth first traversal of the graph.
This string can be parsed and an object with information about the molecule
deduced.

For example, in RDKit, you can use a function like `MolFromSmiles('c1ccccc1')`  
to get a RDKit molecule object. The object contains all kinds of information,
for example you could get the molecular weight, the number of atoms, and much
more. Importantly, you can do structure comparisons, for example checking if
two molecules are equal. RDKit offers a lot of functionality.

## Some background

Database management systems (DBMS) do not usually have chemistry functionality in them.
Extensions like the RDKit Postgres extension allows a user to get cheminformatic
functionality in the DBMS. There are other extensions/cartridges for different
DBMS's offered by different groups.

I've been learning more about DBMS's and duckdb seemed quite interesting to me.
I also like chemistry. I didn't find anything in
a quick google search on a duckdb extension with RDKit, so I thought it would
be interesting to try build one, and hoped that it could be
useful to others who may want to do cheminformatic analyses in duckdb.

Through working on this extension, I've already learned a lot of new things.

I had used RDKit heavily in Python in the past for years, but never dove into
the inner workings. Working on this project brought me to take a closer look
into RDKit's internals on the Molecule object. I also learned more about how
duckdb is working. It's also my first C++ project, which intimidated me.
It took me a while to figure out how to even build duckdb with RDKit using
CMake and make, but eventually I got somewhere and it compiled.

I'm a novice in all these topics, and there's much more to learn.
Nevertheless, I wanted
to write down what I learned so far here.
With that, my disclaimer is that the way I've
implemented things could be/are, the most naive ways to do it.

## Integrating RDKit with duckdb

In order to give duckdb RDKit functionality, I first created a new type for the
RDKit molecule, Mol. Then, I created functions using RDKit to convert from SMILES to Mol
and to serialize the Mol, etc. After that, I wrapped around those functions with
duckdb functions so that that the RDKit functionality can be called in the duckdb
system.

### The Mol type

First, duckdb needs to know what a Molecule is.

The Mol object in RDKit holds a lot of information about the atoms and bonds, etc.
of the molecule. This information is used for further operations on the
molecule. The Mol object for duckdb will be serialized to binary using RDKit's
MolPickler and stored as binary (BLOB) in duckdb.

To create a new type in duckdb, the LogicalType struct is used with the BLOB
enum.

```c++

LogicalType Mol() {
    auto blob_type = LogicalType(LogicalTypeId::BLOB);
    blob_type.SetAlias("Mol");
    return blob_type;
}

```

Then the type needs to be registered

```c++

void RegisterTypes(DatabaseInstance &instance) {
    // Register Mol type
    ExtensionUtil::RegisterType(instance, "Mol", duckdb_rdkit::Mol());
}

```

### Converting between strings to Mol to binary and back

First, I made functions wrapping RDKit, and then called these functions with
code from duckdb that allows this RDKit functionality to work in the duckdb system.

#### Converting from SMILES to Mol

I first tried to implement converting the SMILES to Mol in this way:

```c++
RDKit::ROMol *mol1;
mol1 = RDKit::SmilesToMol("CC");
```

While this worked, when I tried to test the extension on a table from chembl
with about 2.3 million rows with molecular structures, the process running
duckdb would get killed after a while.
I eventually figured out that this way of converting SMILES to Mol causes a memory leak.

I wrote some simple programs and ran them with valgrind to experiment with
how I should properly use RDKit so that I wouldn't get a memory leak.
Running the above with valgrind showed that there were bytes definitely lost
and indirectly lost.
I tried deleting the object after I was done with it, but still ran into different
issues with that approach.

To fix this, I followed the pattern set by chemicalite and used a `std::unique_ptr`
and reset the unique_ptr with the `RDKit::ROMol` pointer that is returned by
`RDKit::SmilesToMol`.

```c++
std::unique_ptr<RDKit::ROMol> mol;
mol.reset(RDKit::SmilesToMol(smiles));
```

After the fact, it seems obvious I should have used smart pointers to begin with,
but I found it helpful to experiment with the approaches and wrap my head around
RDKit and C++.

This removed one source of the memory leak.

#### Serializing and deserialzing Mol

After creating the RDKit Mol, I wanted to serialize the object to binary so
that it could be stored as a BLOB in duckdb. Serializing the object was quite
straightforward, but my first attempt at deserializing it ran into several problems.

One problem was that I had another memory leak when deserialzing the molecule.
The memory leak arose from deserializing the molecule in this way:

```c++

std::string pickle;
RDKit::MolPickler::pickleMol(mol1, pickle);
RDKit::MolPickler::molFromPickle(pickle, mol1);

```

First the molecule object, `mol1` is serialized by `RDKit::MolPickler::pickleMol`
and stored in `pickle`. When I deserialize `pickle` with `molFromPickle`,
I found memory leaks with valgrind. I again tried deleting `mol1` but it didn't
solve my problems. Maybe I did it incorrectly?

The other way to do it, again following chemicalite is as follows:

```c++
  std::unique_ptr<RDKit::ROMol> mol2(new RDKit::ROMol());
  RDKit::MolPickler::molFromPickle(pickle, *mol2);
```

This creates a unique_ptr to a `RDKit::ROMol` allocated on the heap.
Then the mol is deserialized into the new mol object.
However, if I did it like below, I got an empty molecule error.

```c++
std::unique_ptr<RDKit::ROMol> mol;
RDKit::MolPickler::molFromPickle(pickle, *mol2);

```

These changes along with the above solved the memory leak errors.
After building other functions, the extension could process the 2.3 million molecules,
and run some custom functions, which I will describe below, on those molecules
without running out of memory.

### Calling these function from duckdb

The above functions just wrap around RDKit, but are not callable from duckdb yet.
Functions callable by duckdb need to work within the duckdb system. For example,
they may need to know how to accept DataChunks from duckdb and work with duckdb Vectors.

I think of it like putting the RDKit code into a container that duckdb recognizes.

Watching [this talk from Mark Raasveldt on duckdb internals] helped me get a little
better understanding on how duckdb works, but there is a lot I don't quite get yet,
and I'm basically stumbling around to get things working.
I found examples of implementing new functions from other extensions in duckdb,
especially the spatial extension.
Also, people in the duckdb discord have been very helpful in pointing me in the
right direction.

Below is an example of how to make a function that is callable by duckdb.

```c++
void mol_from_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      smiles, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_mol_from_smiles(smiles.GetString());
          auto pickled_mol = rdkit_mol_to_binary_mol(*mol);
          return StringVector::AddString(result, pickled_mol);
        } catch (...) {
          mask.SetInvalid(idx);
          return string_t();
        }
      });
}
```

This will get a SMILES from a DataChunk and then this SMILES as well as a result
Vector is given to a UnaryExecutor. My guess is that this all puts things in a
context or shape that duckdb understands, and plugs into the whole system.

The Executor has a lambda function that takes a `string_t` which is a duckdb
string type, and then here I can call the functions I wrote early that wraps
RDKit functions, `rdkit_mol_from_smiles` to create the molecule object, and then
`rdkit_mol_to_binary_mol` to serialize that object.

#### Bad pickle format problem

I ran into a problem where when I tried to deserialize the binary molecule with
RDKit, I got a bad pickle format error.

Logging the binary that was created from serializing the Mol, I found that every time
I ran the function, I got a different hex representation back.
I struggled with this error for some time, and experimented with simple
RDKit programs without duckdb to check if I was seeing this from RDKit. I also
tested the same SMILES and serialized it in the Postgres RDKit extension.
I did not see this problem there and so it seemed there was something going on
in duckdb that was changing the binary representation.

I asked the people in the duckdb discord and quickly got a response that I need
to return a `string_t` and attach it to a `StringVector` in the UnaryExecutor.
I was returning just `std::string` and not using a `StringVector`. Looking into
`string_t` in duckdb indicated that there are some extra bytes added to the string
by duckdb. It reminded me of the portion of the duckdb internal talk on strings
and how duckdb can quickly filter out strings by attaching some bytes of a string
to the front of this representation, allowing it to quickly check if the first
part of a string matches a query without doing more costly fetching of the
string. (If I understood it correctly)

Anyhow, now with a StringVector being returned, things started working.

### is_exact_match()

It's not interesting to just have a binary molecule in the DB. We want to be able
to do interesting stuff with that binary molecule, for example, finding if
there are certain molecules in our data.

The `is_exact_match(mol1, mol2)` function returns a boolean indicating if the
two molecules being compared are the same.

Long story short, simply comparing strings is not going to work
(see [depth first blog] for more information).
For example, `Cc1ccccc1` and `c1ccccc1C` are the same molecule, toluene.
There's actually several more [issues] to consider when comparing molecules.

To implement is_exact_match, I followed (copied) the chemicalite code to
compare two molecules, which is also similar to the Postgres RDKit extension
code, and I just altered the return type to a boolean.

To create the duckdb function, I have the following:

```c++

static void is_exact_match(DataChunk &args, ExpressionState &state,
                           Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  // args.data[i] is a FLAT_VECTOR
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_blob, string_t &right_blob) {
        // need to first deserialize the blob which is stored in the DB
        std::unique_ptr<RDKit::ROMol> left_mol(new RDKit::ROMol());
        std::unique_ptr<RDKit::ROMol> right_mol(new RDKit::ROMol());

        RDKit::MolPickler::molFromPickle(left_blob.GetString(), *left_mol);
        RDKit::MolPickler::molFromPickle(right_blob.GetString(), *right_mol);
        auto compare_result = mol_cmp(*left_mol, *right_mol);
        return compare_result;
      });

```

This `BinaryExecutor` will take a left and right `string_t`, which should be a
binary molecule, and then deserialize it to a molecule object.

From there, the molecule comparison function, `mol_cmp` can be run to see if the
two molecules are the same.

This works, but in my tests on chembl, it's much slower than the Postgres extension. Adding
an index to the molecule column didn't help, and I confirmed that duckdb is really
fast on a different column, without an index, that is just an integer. This
suggests that there is a problem with my is_exact_match implementation.

For a query `is_exact_match(molecule_column, 'CCC')`, the right side molecule is
not changing.
First, the VARCHAR is converted to a Mol, then to a binary molecule because
the code in BinaryExecutor expects a binary molecule.
Then, I need to deserialize the blob back to a molecule,
and then I run the compare code.
So, for the right side, there would be these steps: VARCHAR -> Mol -> binary -> Mol -> compare,
and this would happen everytime the execute function is ran, which I guess is a
lot if there's 2.3 million tuples.

It would be better to just calculate the molecule once for the unchanging `CCC`,
and use that Mol object for all the comparisons. In a simple program, this is not
difficult; just convert the SMILES outside a for loop and store it in a variable.
Then re-use that.

In duckdb, I'm not sure how to do that.
So I went to the discord and asked for some help. I got a few suggestions of
how I might optimize this code. That is where I left off.

Thanks for reading.

[duckdb_rdkit]: https://github.com/bodowd/duckdb_rdkit
[duckdb]: https://duckdb.org/
[RDKit]: https://www.rdkit.org/docs/Overview.html
[spatial]: https://github.com/duckdb/duckdb_spatial
[chemicalite]: https://github.com/rvianello/chemicalite
[Daylight]: https://www.daylight.com/dayhtml/doc/theory/
[depth first blog]: https://depth-first.com/articles/2019/04/12/the-smiles-substructure-search-fallacy/
[this talk from Mark Raasveldt on duckdb internals]: https://www.youtube.com/watch?v=bZOvAKGkzpQ&ab_channel=CMUDatabaseGroup
[issues]: https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/AM0PR04MB5315B6D8DE8B380CAAF1211FE2AF9@AM0PR04MB5315.eurprd04.prod.outlook.com/
