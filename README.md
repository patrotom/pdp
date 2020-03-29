# MI-PDP.16

This repository contains Parallel and Distributed Programming (MI-PDP.16) course homework from CTU FIT.

* [task1.cpp](src/task1.cpp) – sequential algorithm
* [task2.cpp](src/task2.cpp) – OpenMP task parallelism
* [task3.cpp](src/task3.cpp) – OpenMP data parallelism

## How to compile, execute and run tests

You can take advantage of prepared Makefile recipes, which ease the whole process.

### Compilation

``` bash
make compile${NUMBER}
```

### Compilation for Debugging

``` bash
make debug${NUMBER}
```

### Execute (without input)

Runs compilation automatically.

``` bash
make task${NUMBER}
```

### Run Tests

Runs compilation automatically and uses [tests](tests) directory files as an input.

``` bash
make tests${NUMBER}
```

`${NUMBER}` is a task number you want to run a particular recipe for.

Example (tests for task 2):

``` bash
make tests2
```
