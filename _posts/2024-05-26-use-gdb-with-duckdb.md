---
layout: post
title: use gdb with duckdb
description: ""
summary: ""
tags: [duckdb, databases]
---

I am using GNU 11.3.0. I had some issues with clang

Code is in `examples/embedded-c++/main.cpp`. Will use this to enter the Query function

Built duckdb with debug mode: `GEN=ninja make debug` in the root directory of the project
Go to `examples/embedded-c++` and in CMakeLists.txt use the debug duckdb build:

```git
-link_directories(../../build/release/src)
+link_directories(../../build/debug/src)
```

Then, in the Makefile under `examples/embedded-c++` compile example.cpp with the debug
flag:

```git
 duckdb:
-       cd ../.. && make
+       cd ../.. && GEN=ninja make debug DISABLE_SANITIZER=1

 main:
        mkdir -p build
-       cd build &&  cmake -DCMAKE_BUILD_TYPE=Debug .. && make
+       cd build &&  cmake -DENABLE_SANITIZER=FALSE -DENABLE_UBSAN=0 -DCMAKE_BUILD_TYPE=Debug .. && make
        build/example


```

Now run gdb with `gdb build/example`

(May not need this step anymore after explicitly disabling the sanitizer)
Need to set LD_PRELOAD to libasan shared object in gdb
`(gdb) set environment LD_PRELOAD /usr/lib/x86_64-linux-gnu/libasan.so.6`

I think that if I enable the sanitizer, need to make sure to enable it in both
the compilation of duckdb and in my main.cpp in example. But I didn't try this

Now, set the breakpoint, and run and step into the code.

When I ended up in something like `common/unique_ptr.hpp` or something from the c++ standard library, if I kept stepping through with `s`, it would not go deeper into
the function.
I needed to use `finish` to get out of that function, and then I was able to keep going in the function that I was in
see: https://stackoverflow.com/a/43166815

Example:

- breakpoint at

```c++
auto result = con.Query("SELECT * FROM integers");
```

- run gdb
- I wanted to step into `con.Query` and then inside that function I wanted to continue on
  and step into `context->Query(query, false);`

```gdb
16		auto result = con.Query("SELECT * FROM integers");
(gdb) s
(gdb) n
unique_ptr<MaterializedQueryResult> Connection::Query(const string &query)
91		auto result = context->Query(query, false);
duckdb::shared_ptr<duckdb::ClientContext, true>::operator-> (this=0x7fffffffdc80) at ../duckdb/src/include/duckdb/common/shared_ptr_ipp.hpp:206
206				const auto ptr = internal.get();
(gdb) finish
```

- Here I needed to use `finish` and not `n` when the debugger was in the `shared_ptr_ipp.hpp` code, in order to continue on.
  If I kept using `n` it would not go deeper into the function and into `context->Query`

```gdb
91		auto result = context->Query(query, false);
Value returned is $1 = (duckdb::ClientContext *) 0x61600004bf90
```

- Now I was able to continue stepping deeper into the function

```gdb
(gdb) s
920	unique_ptr<QueryResult> ClientContext::Query(const string &query, bool allow_stream_result) {
```

## extra

### Compiling with clang

I needed to add this to CMakeLists.txt:

```cmake
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")

```

### To get definitions working with clangd in neovim

Run `GEN=ninja make clangd`
