# RAM-MSA

Recursive Anytime Memory-bounded Multiple Sequence Alignment (RAM-MSA) is a CPU-based sequential MSA algorithm.

## Prerequisites

OS: linux / mac OS
c++ boost library

## Compile RAM-MSA

To compile RAM-MSA, first enter the `RAM-MSA` folder. Then,

```
mkdir build
cd build
cmake ..
make
```

## Run RAM-MSA

For the simplest example, enter the `build` folder. Then, 

```
./RAM-MSA
```

This will compute the exact MSA on the `1ac5.fasta` sequence set in the `data` folder with the default settings. 

The following settings are available:

- `-f`: Specific an input FASTA file
- `-t`: Specific the substitution matrix. Default is PAM250.
  + `-t PAM250`: PAM250 substitution matrix (using path cost, int type, gap open penalty = 0, gap extension penalty = -30)
  + (Under construction) `-t BLOSUM62`: BLOSUM62 substitution matrix (using path score, double type, gap open penalty = 1.53, gap extension penalty = 0.00)
- `-m`: Set the memory limit ratio. Default is 0.8, which sets the memory limit to 0.8 * available RAM. Set the ratio to negative to disable the memory-bound strategy.
  + Memory-bound strategy is disabled on mac OS due to some technical issues.

For example, you can compute the exact MSA on `2ack.fasta` without memory-bound strategy using the following command:

```
./RAM-MSA -f "../data/2ack.fasta" -m -1
```

## Other details

1. Output file: `anytime_results.txt`
2. The running time limit is set to 10 million seconds (about 4 months). The search terminates if the timie limit is exceeded.
3. The memory consumption of open & closed nodes are theoretically computed under linear gap penalty conditions. Which means (1) the acutally memory consumption could be different, and (2) it doesn't apply to affine gap penalty cases.
4. The available memory

## Trouble shooting

- CMake Error: CMake 3.15 or higher is required.
  + Solution: In `CMakeLists.txt` line 9, change `VERSION 3.15` to `VERSION 3.14` or lower. Tested with 3.15.
