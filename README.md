## Supplementary Source Code to the Manuscript: Space Net Optimization for Single Objective Bound Constrained Problems.

**This repository holds the supplementary source code to the manuscript entitled "Space Net Optimization for Single Objective Bound Constrained Problems" submitted to a journal for consideration of publication.** 

It contains a C++ implementation of a set of Space Net Optimization (SNO) based algorithms for solving single objective bound-constrained problems. You may refer to a preliminary version of the manuscript (an early draft) to find a brief introduction to the SNO in ArXiv (https://arxiv.org/abs/2306.00043). This repository consists of two programs: [SNO with particle swarm optimization (PSO)](https://github.com/cwtsaiai/sno/blob/main/SNO_PSO) and [SNO with differential evolution (DE)](https://github.com/cwtsaiai/sno/blob/main/SNO_DE), the purpose of which is to make it easier for the audience to understand the basic idea of SNO. Then, based on the transition operators of PSO and DE, the transition operators of SNO are built to show the possibility of SNO. 
## 1. Organization of the programs

- **SNO with PSO**
    - [main_snopso.cpp](https://github.com/cwtsaiai/sno/blob/main/SNO_PSO/main_snopso.cpp): This is where the main function that calls SNO with PSO defined in sno-pso.h for solving single objective bound-constrained problems is defined.
    - [Makefile](https://github.com/cwtsaiai/sno/blob/main/SNO_PSO/Makefile): for compiling and recompiling the program.
    - [sno-pso.h](https://github.com/cwtsaiai/sno/blob/main/SNO_PSO/sno-pso.h): The main source code of SNO with PSO.
    - [test_function.h](https://github.com/cwtsaiai/sno/blob/main/SNO_PSO/test_function.h): The test benchmarks of single objective bound-constrained problems.

- **SNO with DE**
    - [main_snode.cpp](https://github.com/cwtsaiai/sno/blob/main/SNO_DE/main_snode.cpp): This is where the main function that calls SNO with DE defined in sno-de.h for solving single objective bound-constrained problems is defined.
    - [Makefile](https://github.com/cwtsaiai/sno/blob/main/SNO_DE/Makefile): for compiling and recompiling the program.
    - [sno-de.h](https://github.com/cwtsaiai/sno/blob/main/SNO_DE/sno-de.h): The main source code of SNO with DE.
    - [test_function.h](https://github.com/cwtsaiai/sno/blob/main/SNO_DE/test_function.h): The test benchmarks of single objective bound-constrained problems.

## 2. Directory tree
- [src-tree](https://github.com/cwtsaiai/sno/blob/main/src-tree.txt)

## 3. Building and running the programs

1. To clone the repository, type

    ```
    git clone https://github.com/cwtsaiai/sno.git
    ```

2. To compile and/or recompile the program, change to the directory /SNO_PSO or /SNO_DE first and then type

    ```
    make
    ```

3. To run the program, say, SNO with PSO, type either

    ```
    ./main_snopso
    ```

    for the results using default parameters
        
    or
    
    ```
    ./main_snopso [1] [2] [3] [4] [5] [6] [7] [8] [9] [10] [11] [12]
    ```
    
    for the results using user-supplied parameters, such as 

    ```
    ./main_snopso 1 1000 30 10 10 10 10 10 0.2 0.2 10 1
    ```
    
    The meanings of the parameters are shown below.

    [1] The number of rounds (runs) to carry out the test function.

    [2] The maximum number of evaluations each run.

    [3] The number of dimensions of the test function.

    [4] The size of the space net.

    [5] The initial size of the explorer.

    [6] The final size of the explorer.

    [7] The initial size of the miner.

    [8] The final size of the miner.

    [9] The control parameter $F$ for SNO with DE. The control parameter $w$ for SNO with PSO.

    [10] The control parameter $Cr$ for SNO with DE. The control parameter $c_1$ for SNO with PSO.

    [11] The maximum number of space net adjustments for each search history. 

    [12] Index of the test function, as specified in test_function.h.

## 4. LICENSE

**Copyright (c) 2022-2024 Chun-Wei Tsai and AILAB@NSYSU**  
**All rights reserved.**  

**This work (the supplementary code in this repository) is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](http://creativecommons.org/licenses/by-nc-nd/4.0/).**
