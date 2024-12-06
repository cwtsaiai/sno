## Supplementary Source Code to the Manuscript: Space Net Optimization for Single Objective Bound Constrained Problems.

**This repository holds the supplementary source code to the manuscript entitled "Space Net Optimization for Single Objective Bound Constrained Problems," which has now been submitted to a journal.** 

It contains an implementation of a set of Space Net Optimization (SNO) based algorithms in C++ for solving single objective bound-constrained problems. You might be able to refer to a preliminary version of the manuscript (early draft) to find a brief introduction to the SNO in ArXiv (https://arxiv.org/abs/2306.00043). This repository can be divided into twofold: [SNO with particle swarm optimization (PSO)](https://github.com/cwtsaiai/sno/blob/SNO_PSO) and [SNO with differential evolution (DE)](https://github.com/cwtsaiai/sno/SNO_DE). It is used to let the audience easily understand the basic idea of SNO, and then, based on the transition operators of PSO and DE, the transition operators of SNO are designed to show the possibility of SNO. 
## 1. Organization of the programs
- **SNO with PSO**
 - [main_snopso.cpp] 
(https://github.com/cwtsaiai/sno/blob/SNO_PSO/main_snopso.cpp): The main function file to call the sno-pso.h to run SNO with PSO for  solving single objective bound-constrained problems
 - [makefile](https://github.com/cwtsaiai/sno/blob/SNO_PSO/makefile): The file that can be use to compile the source code.
 - [sno-pso.h] (https://github.com/cwtsaiai/sno/blob/SNO_PSO/sno-pso.h): The main source code of SNO with PSO.
 - [test_function.h] (https://github.com/cwtsaiai/sno/blob/SNO_PSO/test_function.h): The test benchmarks of single objective bound-constrained problems.
  
- **SNO with DE**
 - [main_snode.cpp] 
(https://github.com/cwtsaiai/sno/blob/SNO_DE/main_snode.cpp): The main function file to call the sno-de.h to run SNO with DE for  solving single objective bound-constrained problems
 - [makefile](https://github.com/cwtsaiai/sno/blob/SNO_DE/makefile): The file that can be use to compile the source code.
 - [sno-de.h] (https://github.com/cwtsaiai/sno/blob/SNO_DE/sno-de.h): The main source code of SNO with DE.
 - [test_function.h] (https://github.com/cwtsaiai/sno/blob/SNO_DE/test_function.h): The test benchmarks of single objective bound-constrained problems.


## 2. Directory tree
- [src-tree](https://github.com/cwtsaiai/sno/src-tree.txt)

## 3. Building and running the programs
1. To clone the repository, type
```
git clone https://github.com/cwtsaiai/sno.git
```

2. To compile the C++ code, change to the directory /SNO_PSO or /SNO_DE first and then type
```
make
```

3. To run a program, say, SNO with PSO, type either
```
./main_snopso
```
    to get the results with the default settings,

    or
```
./main_snopso [1] [2] [3] [4] [5] [6] [7] [8] [9] [10] [11] [12]
```
    to get the results with the predefined settings by the user, which means we use preferred parameter settings for specific test functions, such as 
```
./main_sno 1 1000 30 10 10 10 10 10 0.2 0.2 10 1
```
The meanings of these parameters are shown below.

[1] The number of rounds (runs) for the test function will be carried out.

[2] The number of maximum evaluations each run.

[3] The number of dimensions of the test function.

[4] The size of the space net.

[5] The initial size of the explorer.

[6] The final size of the explorer.

[7] The initial size of the miner.

[8] The final size of the miner.

[9] The control parameter $F$ in DE version. The control parameter $\omega$ in PSO version.

[10] The control parameter $Cr$ in DE version. The control parameter $c_1$ in PSO version.

[11] The maximum number of space net adjustments for each search history. 

[12] Specific index of test function. Test_function.h contains ten different functions.

## 4. LICENSE

**Copyright (c) 2022-2023 Chun-Wei Tsai and AILAB@NSYSU**  
**All rights reserved.**  

**This work (the supplementary code in this repository) is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](http://creativecommons.org/licenses/by-nc-nd/4.0/).**
