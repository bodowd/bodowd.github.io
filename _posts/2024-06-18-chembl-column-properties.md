---
layout: post
title: chembl column properties analysis
description: ""
summary: ""
tags: [chembl, chemistry, rdkit, databases]
---

In `An Empirical Evaluation of Columnar Storage Formats` by X. Zeng, et al., VLDB 2023 [[1]],
the authors studied the Parquet and ORC file formats in order to inform the design of
the next generation of columnar storage formats. They began by analyzing the
parameter distribution of data sets, and used the insights to construct a benchmark
against which they could test the file formats in a variety of tasks.

After reading the paper, I wondered what the parameter distribution of a
chemical life-science database would look like, and how it would compare to the
distributions the authors observed.

The [chembl] database contains 2.4 million compounds and all
kinds of data associated with the compounds, like assay results and drug targets.

I don't know how representative the chembl database is of the databases in chemical life science
companies (pharma, agriculture), but perhaps it is similar enough.

In the paper, the authors calculate the following properties on int, float, and
string type columns from a collection of datasets: number of distinct values ratio (NDV),
null ratio, value range, sortedness, and skew pattern.

I looked at NDV, null ratio, and value range properties.

I also took a look at chemistry specific data, like molecule objects from RDKit,
which is necessary to perform cheminformatic workloads.

Jump to the conclusion [here](#conclusion).
Code is [here]

# Summary of results

## NDV

NDV is calculated by the number of distinct values in the column divided by the number
of values in the column.

### Int

| ndv            | % of dataset | cumulative |
| :------------- | -----------: | ---------: |
| (0.0, 1e-05]   |     0.114286 |   0.114286 |
| (1e-05, 1e-04] |    0.0571429 |   0.171429 |
| (1e-04, 1e-03] |    0.0619048 |   0.233333 |
| (1e-03, 1e-02] |    0.0714286 |   0.304762 |
| (1e-02, 1e-01] |     0.104762 |   0.409524 |
| (1e-01, 1]     |     0.590476 |          1 |

In chembl ~30% of columns have a NDV ratio <= 0.01 while in the paper's analysis
~80% of integer columns had an NDV ratio < 0.01.
Since my analysis looks directly at the database tables' columns, I have a lot
of primary keys and foreign keys which are just increasing sequences.
These columns would have high NDV ratios.
I don't know what the datasets from the paper look like, but if these datasets
are flattened or denormalized, they could have less `id` columns from primary keys or foreign keys.

### Float

| ndv            | % of dataset | cumulative |
| :------------- | -----------: | ---------: |
| (0.0, 1e-05]   |            0 |          0 |
| (1e-05, 1e-04] |    0.0322581 |  0.0322581 |
| (1e-04, 1e-03] |      0.16129 |   0.193548 |
| (1e-03, 1e-02] |     0.322581 |   0.516129 |
| (1e-02, 1e-01] |     0.354839 |   0.870968 |
| (1e-01, 1]     |    0.0967742 |   0.967742 |

In the paper, the authors found 60% of floating-point columns had a NDV ratio
< 0.1, which they classified as significant value repetitions.
In the chembl database, I found ~87% of floating-point columns had a NDV ratio
<= 0.1, even higher than the datasets the authors looked at.
As they point out in the paper, Dictionary Encoding would be effective here.

I suspect that this result could be due to values in the database that come
from assays that are more discrete than one would expect, but I didn't investigate
this further in the chembl dataset since this could be assay dependent.

### String

| ndv            | % of dataset | cumulative |
| :------------- | -----------: | ---------: |
| (0.0, 1e-05]   |    0.0721649 |  0.0721649 |
| (1e-05, 1e-04] |    0.0618557 |   0.134021 |
| (1e-04, 1e-03] |    0.0515464 |   0.185567 |
| (1e-03, 1e-02] |     0.113402 |   0.298969 |
| (1e-02, 1e-01] |      0.14433 |   0.443299 |
| (1e-01, 1]     |     0.546392 |   0.989691 |

In the chembl dataset, only ~30% of the string columns have a NDV ratio <=0.01,
as opposed to ~60% of the string columns in the datasets in the paper.

## Null Ratio

### Int

| null_ratio     | % of dataset | cumulative |
| :------------- | -----------: | ---------: |
| 0              |     0.809524 |   0.809524 |
| (0.0, 1e-05]   |            0 |   0.809524 |
| (1e-05, 1e-04] |            0 |   0.809524 |
| (1e-04, 1e-03] |   0.00952381 |   0.819048 |
| (1e-03, 1e-02] |   0.00952381 |   0.828571 |
| (1e-02, 1e-01] |    0.0857143 |   0.914286 |
| (1e-01, 1]     |    0.0857143 |          1 |

~91% of the int columns have a null ratio <= 0.1 (not many nulls).
This seems similar to the findings in the paper.

### Float

| null_ratio     | % of dataset | cumulative |
| :------------- | -----------: | ---------: |
| 0              |    0.0967742 |  0.0967742 |
| (0.0, 1e-05]   |    0.0645161 |    0.16129 |
| (1e-05, 1e-04] |            0 |    0.16129 |
| (1e-04, 1e-03] |            0 |    0.16129 |
| (1e-03, 1e-02] |    0.0322581 |   0.193548 |
| (1e-02, 1e-01] |     0.419355 |   0.612903 |
| (1e-01, 1]     |     0.387097 |          1 |

\~61% of the float columns have a null ratio <= 0.1, and \~39% of columns have a
null ratio > 0.1 , which seems lower than the findings in the paper. They
found ~85% of the float columns have a null ratio < 0.1.

### String

| null_ratio     | % of dataset | cumulative |
| :------------- | -----------: | ---------: |
| 0              |     0.570447 |   0.570447 |
| (0.0, 1e-05]   |   0.00343643 |   0.573883 |
| (1e-05, 1e-04] |   0.00343643 |    0.57732 |
| (1e-04, 1e-03] |     0.024055 |   0.601375 |
| (1e-03, 1e-02] |      0.04811 |   0.649485 |
| (1e-02, 1e-01] |    0.0652921 |   0.714777 |
| (1e-01, 1]     |     0.285223 |          1 |

~71% of the string columns have a null ratio <= 0.1. Unlike the findings in the paper,
the chembl database has less nulls in the string columns than the float columns,
whereas the authors found that strings had more null values than the other types.

## Value Range

### Int

| range      | % of dataset | cumulative |
| :--------- | -----------: | ---------: |
| [0.0, 1]   |     0.119048 |   0.119048 |
| (1, 1e1]   |    0.0857143 |   0.204762 |
| (1e1, 1e2] |    0.0904762 |   0.295238 |
| (1e2, 1e3] |    0.0857143 |   0.380952 |
| (1e3, 1e4] |     0.109524 |   0.490476 |
| (1e4, 1e5] |     0.147619 |   0.638095 |
| (1e5, inf] |     0.361905 |          1 |

~64% of int columns have a value range <= 100,000 and ~36% > 100,000

### String

| range       | % of dataset | cumulative |
| :---------- | -----------: | ---------: |
| [0.0, 5]    |     0.340278 |   0.340278 |
| (5, 10]     |    0.0763889 |   0.416667 |
| (10, 25]    |     0.138889 |   0.555556 |
| (25, 50]    |      0.09375 |   0.649306 |
| (50, 100]   |     0.142361 |   0.791667 |
| (100, 1000] |     0.149306 |   0.940972 |
| (1000, inf) |    0.0590278 |          1 |

The length of a string is measured in bytes. ~94% of the string data is less than
or equal to 1,000 bytes.

### In the context of cheminformatics workloads

I took a deeper look at the length of strings and found there are some big
string records. This include description of assays, but also there are molecular file
formats in this database, which is quite common to chemistry data.

I found a molfile can reach a length of 65,000 bytes.

In the context of this paper and study, the focus is on looking at the data from a perspective of
the compression and fast decompression of the data for use in analytical database systems.
In cheminformatics workloads, we would not typically work on the raw string.

At least in RDKit, you typically want the molecule object, a particular encoding of
the information about a molecule, which is constructed by parsing a molecule
format, like the molfile or a SMILES string. With this molecule object, cheminformatic
operations can be performed. For example, if you want to find if a molecule, or
similar molecule, exists in the database, typically you would run the algorithm
comparing data between molecule objects and not between string representations.

Parsing the string format into a molecule object, is
the most computationally expensive part, at least in my experiments. Then we need
to serialize the object to binary for storage. This serialization and later deserialization
seems to be cheaper compared to first part of parsing the string into the molecule object.
If computation needs to be performed on these molecule objects, it would make
sense to pre-compute the molecule objects and store these binary encodings
in the database. I thought that the size of these binary encodings might be
of interest in the context of this paper and analysis.

I created a molecule object column with the RDKit Postgres extension
in the `compound_structures` table in the database which contains chemical structures.
I took a look at the size of these objects and found the max byte length for the
molecules in that table was 15,762, the min byte length was 65, and the average was
454.9 bytes.

The max canonical SMILES (string representation) byte length was 2045 bytes,
the min was 1, and the average was 57, indicating that the molecule object is a lot bigger than the strings.
That seems to make sense to me since there is a lot of information about the molecule's
structure and properties in the molecule object.

# <a name="conclusion"></a>Conclusion

Although some of the column properties look different in the chembl database compared
to the datasets analyzed in the paper, I'm not sure if it would substantially
change the insights regarding the performance of Parquet and ORC files.

I think the biggest difference is regarding the molecule information,
namely the molecule objects needed for cheminformatic operations.
I am curious how an analytical database system could leverage the
binary molecule objects in its query execution. I am also curious about how
these binary objects could affect the performance of Parquet and ORC files,
and future formats, if at all.

[1]: https://15721.courses.cs.cmu.edu/spring2024/papers/02-data1/p148-zeng.pdf
[chembl]: https://www.ebi.ac.uk/chembl/
[here]: https://github.com/bodowd/chembl-column-analysis
