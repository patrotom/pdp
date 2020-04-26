# Analysis

This is an analysis of the execution of particular tasks.

* Sequential – [task1.cpp](../src/task1.cpp)
* Task parallelism – [task2.cpp](../src/task2.cpp)
* Data parallelism – [task3.cpp](../src/task3.cpp)
* MPI – [task4.cpp](../src/task4.cpp)

## Execution Times

### Local

Times are calculated using `std::chrono` library. All times are represented in seconds.

* [task4.cpp](../src/task4.cpp) was run on 4 processes.
* Configuration:
  * Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz
  * 8 CPUs

| Algorithm\Input      | mvr_20_10_5 | mvr_20_15_8 | mvr_30_10_10 | mvr_35_20_10 | mvr_40_20_10 | mvr_45_20_15 | mvr_45_25_15 | mvr_45_25_18 | mvr_48_25_18 | mvr_50_20_15 | mvr_51_20_15 |
|----------------------|-------------|-------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|
| **Sequential**       | 0.003       | 0.0009      | 0.041        | 1.61         | 11.567       | 58.066       | 54.868       | 16.518       | 18.043       | 179.601      | 2108.55      |
| **Task parallelism** | 0.018       | 0.012       | 0.04         | 0.825        | 8.015        | 33.969       | 28.391       | 9.176        | 10.452       | 115.541      | 1310.37      |
| **Data parallelism** | 0.006       | 0.006       | 0.03         | 0.979        | 11.616       | 26.334       | 23.923       | 6.832        | 9.076        | 107.823      | 1208.04      |
| **MPI**              | 0.904       | 0.92        | 0.686        | 1.808        | 13.867       | 28.027       | 27.81        | 8.381        | 9.298        | 128.702      | 1718.34      |

### Star

| Algorithm\Input      | mvr_20_10_5 | mvr_20_15_8 | mvr_30_10_10 | mvr_35_20_10 | mvr_40_20_10 | mvr_45_20_15 | mvr_45_25_15 | mvr_45_25_18 | mvr_48_25_18 | mvr_50_20_15 | mvr_51_20_15 |
|----------------------|-------------|-------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|
| **Sequential**       |        |       |         |          |        |        |        |        |        |       |       |
