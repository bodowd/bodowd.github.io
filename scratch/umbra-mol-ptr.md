- Added some print statements to see if this function gets called when
  running operations of VARCHAR data when running duckdb with a file attached
- Idea: since it's all just binary data in the end, convert to `string_t` so that
  duckdb can handle pointer swizzling. When bringing the data into memory, convert
  from `string_t` to an UmbraMol, `umbra_mol_t`. The main difference between the
  two types is the size of the prefix. Then, set the `umbra_mol_t` pointer to the
  `string_t` pointer to the data.

Store the data in a pointer. Convert to and from string_t and let duckdb
handle the pointer swizzling.

With just count prefix in the string_t.

Exact match queries

<!-- | Query | Standard method (s) | Umbra-mol 20-byte prefix (s) | Umbra-mol 10-byte prefix (s) | Umbra-mol as string_t (s) | Postgres control (s) | -->
<!-- | :---- | :------------------ | :--------------------------- | :--------------------------- | :------------------------ | :------------------- | -->
<!-- | 1     | 17.238              | 0.496                        | 0.311                        | 0.141                     | 0.084                | -->
<!-- | 2     | 12.555              | 0.473                        | 0.273                        | 0.103                     | 233                  | -->
<!-- | 3     | 22.196              | 0.364                        | 0.592                        | 0.240                     | 0.162                | -->
<!-- | 4     | 12.245              | 0.359                        | 0.350                        | 0.217                     | 0.900                | -->

| Query | Standard method (s) | Umbra-mol 20-byte prefix (s) | Umbra-mol as string_t (s) | Speedup (string_t vs standard method) |
| :---- | :------------------ | :--------------------------- | :------------------------ | :------------------------------------ |
| 1     | 17.238              | 0.496                        | 0.141                     | 122x                                  |
| 2     | 12.555              | 0.473                        | 0.103                     | 122x                                  |
| 3     | 22.196              | 0.364                        | 0.240                     | 92x                                   |
| 4     | 12.245              | 0.359                        | 0.217                     | 56x                                   |

With counts + dalke fingerprints (12 byte prefix)

Exact matches:

| Query | Standard method (s) | Umbra-mol 20-byte prefix (s) | Umbra-mol (12 byte count + dalke fp) (s) | Postgres control (s) |
| :---- | :------------------ | :--------------------------- | :--------------------------------------- | :------------------- |
| 1     | 17.238              | 0.496                        | 0.211                                    | 0.084                |
| 2     | 12.555              | 0.473                        | 0.175                                    | 233                  |
| 3     | 22.196              | 0.364                        | 0.699                                    | 0.162                |
| 4     | 12.245              | 0.359                        | 0.351                                    | 0.900                |

Substructure matches:

| Query | Standard method (s) | Umbra-mol 10-byte prefix + dalke fp (s) | Umbra-mol (12-byte count + dalke fp) (s) | Postgres (s) |
| :---- | :------------------ | :-------------------------------------- | :--------------------------------------- | ------------ |
| 1     | 14.487              | 19.95                                   | 8.904                                    | 130.59       |
| 2     | 13.96               | 19.47                                   | 0.768                                    | 1.6          |

```sql
SELECT umbra_mol from molecule where umbra_is_substruct(umbra_mol, 'O=CNCCc1ccccc1');
```

133531 molecules. 9.061 seconds

```sql
SELECT rdkit_mol from molecule where is_substruct(rdkit_mol, 'O=CNCCc1ccccc1');
```

133531 molecules. 16.374 seconds

postgres

```sql
select rdkit_mol from compound_structures where rdkit_mol@>'O=CNCCc1ccccc1';
```

133531 molecules. 143 seconds.
