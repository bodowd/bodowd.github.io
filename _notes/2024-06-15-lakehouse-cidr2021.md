---
layout: post
title: Lakehouse (M. Armbrust, et al., CIDR 2021)
description: ""
summary: ""
tags: [databases, notes, lakehouse]
---

Lakehouse: A New Generation of Open Platforms that Unify Data Warehousing and Advanced Analytics (M. Armbrust, et al., CIDR 2021)

### Overview:

- The first generation data analysis platforms were data warehouses which had some problems: coupled compute and storage on premises, and they did not handle unstructured data like video, audio, and text documents
- The second generation platforms used data lakes which provide low cost direct access storage, but still required data warehouses to provide DBMS capabilities, resulting in complex two-tiered architectures with four major problems:
  1. Reliability: difficult to keep the two systems consistent and each ETL step brings additional risks and failure modes
  2. Data staleness: data in the warehouse is stale compared to the data in the data lake
  3. Limited support for advanced analytics: ML systems need to process large datasets with complex non-SQL code
  4. Total cost of ownership: expensive and hard to migrate to different data warehouse solutions because of proprietary formats
- The Lakehouse architecture proposed in this paper combines the advantages of the data lake (low cost direct access storage) and data warehouse (traditional analytical DBMS features like ACID transactions, query optimization, etc…), and a system was built and evaluated to test if this design could be an effective solution.

### Key takeaways:

- The key pieces of the Lakehouse architecture are:
  1. object store like S3
  2. open file format like Apache Parquet
  3. Transactional metadata layer on top of the object store
  4. additional optimizations to improve SQL performance like caching, auxiliary data structures like indexes and statistics, and data layout optimizations
- They report comparable or better performance as the data warehouses they tested against, and also were able to achieve this at a lower cost

### System evaluated:

- Delta Engine: a C++ execution engine for Apache Spark built at Databricks

### Workloads/Benchmarks

- They tested the test duration and cost to run the tests in TPC-DS at a scale factor of 30,000 with four widely used cloud data warehouses on AWS, Azure, and GCP
