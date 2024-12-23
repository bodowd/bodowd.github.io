---
layout: post
title: DuckDB Pragma statements and Client Context note
description: ""
summary: ""
tags: [duckdb]
---

This paper [The Key to Effective UDF Optimization:
Before Inlining, First Perform Outlining](https://www.vldb.org/pvldb/vol18/p1-arch.pdf)
implements its ideas in DuckDB as an extension and I was curious
specifically on how the authors access the SQL UDF to start working on it.

I got the code from Sam Arch's [github](https://github.com/SamArch27/PRISM/tree/main), and
I followed the instructions there on building and running their extension.

In `udf_transpiler_extension.cpp` in the `LoadInternal` function, the
`PragmaFunction::PragmaCall` for `UdfTranspilerPragmaFun` is registered.

```c++
inline String UdfTranspilerPragmaFun(ClientContext &context,
                                     const FunctionParameters &parameters) {
  std::cout << "----- Running UdfTranspilerPragmaFun -----" << std::endl;
  auto udfString = parameters.values[0].GetValue<String>();

  return CompilerRun(udfString);
}

auto udf_transpiler_pragma_function = PragmaFunction::PragmaCall(
  "transpile", UdfTranspilerPragmaFun, {LogicalType::VARCHAR});
ExtensionUtil::RegisterFunction(instance, udf_transpiler_pragma_function);


```

Then, to see when this was called, I traced back to the `ClientContext`.
The context will first call `Query` and then that will call `ParseStatements`.
When the statement is parsed by `ParseStatementsInternal`, it will call
`HandlePragmaStatements`.

```cpp
void PragmaHandler::HandlePragmaStatements(ClientContextLock &lock, vector<unique_ptr<SQLStatement>> &statements) {
	// first check if there are any pragma statements
	bool found_pragma = false;
	for (idx_t i = 0; i < statements.size(); i++) {
		if (statements[i]->type == StatementType::PRAGMA_STATEMENT ||
		    statements[i]->type == StatementType::MULTI_STATEMENT) {

			std::cout << "statements[i]->type: " << StatementTypeToString(statements[i]->type) << std::endl;
			std::cout << "statements query: " << statements[i]->query << std::endl;
			found_pragma = true;
			break;
		}
	}

	std::cout << "!!!!!!!! Found pragma: " << found_pragma << std::endl;
	if (!found_pragma) {
		// no pragmas: skip this step
		return;
	}
	context.RunFunctionInTransactionInternal(lock, [&]() { HandlePragmaStatementsInternal(statements); });
}

```

The UDF is detected by duckdb as a PRAGMA statement.

Here we see the statement has the `transpile` function is passed the UDF.

```shell
statements[i]->type: PRAGMA
statements query: pragma transpile('CREATE FUNCTION addAbs(val1 INT, val2 INT) RETURNS INT AS $$ BEGIN IF val1 < 0 THEN val1 = -val1; END IF; IF val2 < 0 THEN val2 = -val2; END IF; RETURN val1 + val2; END; $$ LANGUAGE PLPGSQL;');
!!!!!!!! Found pragma: 1
----- Running UdfTranspilerPragmaFun -----

```
