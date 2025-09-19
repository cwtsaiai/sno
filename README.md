## A Simplified Version of Space Net Optimization for Single Objective Bound Constrained Problems.

**This repository holds the supplementary source code to the manuscript entitled "Space Net Optimization for Single Objective Bound Constrained Problems" published in the IEEE Transactions on Evolutionary Computation.** 

It contains a C++ implementation of a set of simplified version of Space Net Optimization (SNO) based algorithms for solving single objective bound-constrained problems. You may refer to a preliminary version of the manuscript (an early draft) to find a brief introduction to the SNO in ArXiv (https://arxiv.org/abs/2306.00043), as well as a full version in the IEEE Transactions on Evolutionary Computation (TEVC) (https://ieeexplore.ieee.org/abstract/document/11126067).

This repository consists of two programs: [SNO with differential evolution (DE)](https://github.com/cwtsaiai/sno/tree/main/src/SNO_DE) and [SNO with particle swarm optimization (PSO)](https://github.com/cwtsaiai/sno/tree/main/src/SNO_PSO), the purpose of which is to make it easier for the audience to understand the basic idea of SNO. In these implementations, we adapt the transition operators of DE (SNO-DE) and PSO (SNO-PSO) for the space net to generate portion of the new candidate solutions during the convergence process, thereby demonstrating the impact of the space net, when compared with conventional DE and PSO. The detailed results and discussions can also be seen in Section III.B, Preliminary Tests of the manuscript of TEVC. You may also find [Tables](https://github.com/cwtsaiai/sno/tree/main/notation.pdf) showing the mapping between notation in the paper and variables in the code as well as the mapping between SNO operators in the paper and Functions in the code. Note that the population for the results of the Preliminary Tests described in Table II is of size 120, but SNO-DE and SNO-PSO can also obtain good performance for the population of size 190, too.

The main contributors to these implementations are: Chun-Wei Tsai (who presents the basic idea of SNO, builds the very first version of the prototype to prove the concept, and participates in the tasks of optimizing and organizing the code in the last stage), Yi-Cheng Yang (who builds the early version of the code), Wei-Hong Wang (who fine-tunes the code to speed it up in the later stage), Kung-Ming Liu (who attempts to visualize the results in high-dimensional space in the later stage), Chun-Chih Chao (who participates in the tasks of verifying the code in the last stage), and Ming-Chao Chiang (who reviews, reorganizes, and optimizes all the code in the last stage).

## 1. Organization of the programs

- **SNO with DE**
    - [main_snode.cpp](https://github.com/cwtsaiai/sno/blob/main/src/SNO_DE/main_snode.cpp): This is where the main function that calls SNO with DE defined in sno-de.h for solving single objective bound-constrained problems is defined.
    - [Makefile](https://github.com/cwtsaiai/sno/blob/main/src/SNO_DE/Makefile): for compiling and recompiling the program.
    - [sno-de.h](https://github.com/cwtsaiai/sno/blob/main/src/SNO_DE/sno-de.h): The main source code of SNO with DE.
    - [test_function.h](https://github.com/cwtsaiai/sno/blob/main/src/SNO_DE/test_function.h): The test benchmarks of single objective bound-constrained problems.
    
- **SNO with PSO**
    - [main_snopso.cpp](https://github.com/cwtsaiai/sno/blob/main/src/SNO_PSO/main_snopso.cpp): This is where the main function that calls SNO with PSO defined in sno-pso.h for solving single objective bound-constrained problems is defined.
    - [Makefile](https://github.com/cwtsaiai/sno/blob/main/src/SNO_PSO/Makefile): for compiling and recompiling the program.
    - [sno-pso.h](https://github.com/cwtsaiai/sno/blob/main/src/SNO_PSO/sno-pso.h): The main source code of SNO with PSO.
    - [test_function.h](https://github.com/cwtsaiai/sno/blob/main/src/SNO_PSO/test_function.h): The test benchmarks of single objective bound-constrained problems. 

## 2. Directory tree

- [src-tree](https://github.com/cwtsaiai/sno/blob/main/src-tree.txt)

## 3. Building and running the programs

1. To clone the repository, type

    ```
    git clone https://github.com/cwtsaiai/sno.git
    ```

2. To compile and/or recompile the program, change to the directory sno/src/SNO_PSO or sno/src/SNO_DE first and then type

    ```
    make
    ```

3. To run the program, say, SNO with PSO, type either

    ```
    ./main_snode
    ```

    for the results using default parameters
        
    or
    
    ```
    ./main_snode [1] [2] [3] [4] [5] [6] [7] [8] [9] [10] [11] [12]
    ```
    
    for the results using user-supplied parameters, such as 

    ```
    ./main_snode 1 1000 30 10 10 10 10 10 0.2 0.2 10 1
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

## 4. Visualization of SNO for Some Benchmarks

Here, animations are given to show the isometric and top views of SNO in solving the benchmarks (namely, Ackley, Bent Cigar, Griewank, Rastrigin, and Rosenbrock), thus making it easy to see how the space net is adjusted to fit the landscape of the solution space during the convergence process. 

<div align="left">
<img src="https://github.com/cwtsaiai/sno/blob/main/gif/Ackley.gif" width="400" />
<img src="https://github.com/cwtsaiai/sno/blob/main/gif/BentCigar.gif" width="400"/><br>
<img src="https://github.com/cwtsaiai/sno/blob/main/gif/Griewank.gif" width="400" />
<img src="https://github.com/cwtsaiai/sno/blob/main/gif/Rastrigin.gif" width="400" /><br>
<img src="https://github.com/cwtsaiai/sno/blob/main/gif/Rosenbrock.gif" width="400" />
</div>

## 5. Citation

[Chun-Wei Tsai](https://sites.google.com/site/cwtsai0807/chun-wei-tsai),[Yi-Cheng Yang](), [Tzu-Chieh Tang](), [Che-Wei Hsu](), ["Space Net Optimization for Single Objective Bound Constrained Problems"](https://www.doi.org/10.1109/TEVC.2025.3599438), <i>IEEE Transactions on Evolutionary Computation (TEVC)</i>, 2025, In Press.

## 6. LICENSE

**Copyright (c) 2022-2024 Chun-Wei Tsai and AILAB@NSYSU** 
**All rights reserved.** 

**This work (the supplementary code in this repository) is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](http://creativecommons.org/licenses/by-nc-nd/4.0/).**
