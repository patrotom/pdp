# MI-PDP.16

This repository contains Parallel and Distributed Programming (MI-PDP.16) course homework from CTU FIT.

* [task1.cpp](src/task1.cpp) – sequential algorithm
* [task2.cpp](src/task2.cpp) – OpenMP task parallelism
* [task3.cpp](src/task3.cpp) – OpenMP data parallelism
* [task4.cpp](src/task4.cpp) – MPI

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

### Notes

* There are prepared only `compile` and `tests` recipes for [task4.cpp](src/task4.cpp).
* By default, [task4.cpp](src/task4.cpp) is executed using `mpic++` utility.
  * Default number of MPI processes (-np) is `3`.
  * To change this value, run the makefile `tests` recipe like this:

      ``` bash
      make NP=<number_of_processes> tests4
      ```

  * For example, to use 4 processes you run it like this:

      ``` bash
      make NP=4 tests4
      ```

* By default, [task4.cpp](src/task4.cpp) is executed without defining `OMP_NUM_THREADS`.
  * To change this behavior, follow the steps from the previous example, e.g.:

    ``` bash
    make OMP_NUM_THREADS=<number_of_threads> tests4
    ```

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
