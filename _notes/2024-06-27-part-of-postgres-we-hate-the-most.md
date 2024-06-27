---
layout: post
title: The part of PostgreSQL we hate the most (Pavlo & Zhang)
description: ""
summary: ""
tags: [notes, databases, postgres, mvcc]
---

The Part of PostgreSQL We Hate the Most by Pavlo and Zhang [[1]]

### Overview:

- MVCC is a transaction management technique that enables multiple queries to read
  older versions of rows without getting blocked by transactions doing writes
- Postgres's MVCC implementation has some important consequences on query performance,
  memory usage, and disk space.
- Four problems due to Postgres's implementation of MVCC are discussed
  1. Version Copying
  2. Table Bloat
  3. Secondary Index Maintenance
  4. Vacuum Management

### Key takeaways:

- There are three important design decisions for implementing MVCC:

  1. How to store updates to existing rows
  2. How to find the correct version of a row at runtime
  3. How to remove expired versions that are no longer visible to transactions

- Postgres's design causes the following issues which appear especially on
  write heavy workloads:
  - Version Copying
    - Postgres updates existing rows by copying the entire tuple being modified,
      including unmodified columns to a new slot in a page, could be a different page
      which causes Postgres to have higher memory and disk usage than other systems,
      translating to slower queries and higher cloud costs
  - Table bloat
    - Copying the entire tuple causes table bloat
    - Uses more storage space and can cause dead tuples to accumulate faster
      than garbage collection (vacuum) can clean them up if there are write heavy
      workloads
    - Query performance slows because this requires more disk IO
  - Secondary index maintenance
    - All indexes in the table that the modified tuple belongs to are updated and
      indexes point to physical addresses, as opposed to logical identifiers like
      tuple ids or primary key (faster writes but makes secondary index reads slower)
    - More indexes to update, slower writes, increased cloud costs
  - Vacuum management
    - blocked by long running transactions which results in more dead tuples and stale
      stats which cause a cycle of more long running transactions

[1]: https://www.cs.cmu.edu/~pavlo/blog/2023/04/the-part-of-postgresql-we-hate-the-most.html
