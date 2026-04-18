# RAM-MSA

[![DOI](https://zenodo.org/badge/1096295403.svg)](https://doi.org/10.5281/zenodo.17607120)

Recursive Anytime Memory-bounded Multiple Sequence Alignment (RAM-MSA) is a CPU-based sequential MSA algorithm.

## Prerequisites

- OS: linux / macOS
  + Tested on Ubuntu 22.04.5 LTS, macOS 12.7.6, and Windows 11 x86_64
- CMake (https://cmake.org/download/)
- c++ compiler
- c++ boost library (already prepared in include/boost)

## Compile RAM-MSA

To compile RAM-MSA, first enter the `RAM-MSA` folder. Next,

```
mkdir build
cd build
cmake ..
```

Then

`make` or `ninja`, depending on your environment


## Run RAM-MSA

For the simplest example, enter the `build` folder. Then, 

```
./RAM-MSA
```

This will compute the exact MSA on the SABRE `twi_009` sequence set in the `data` folder with the default settings. 

The following settings are available:

- `-f`: Specify an input FASTA file
- `-rf`: Specify an reference FASTA file. The program will output the sum-of-pairs (SP) accuracy and total column (TC) accuracy if the reference file is valid.
- `-t`: Specify the substitution matrix. Default is PAM250.
  + `-t PAM250`: PAM250 substitution matrix (using path cost, int type, gap open penalty = 0, gap extension penalty = 12)
  + `-t BLOSUM62`: BLOSUM62 substitution matrix (using path score, float type, gap open penalty = 9.5, gap extension penalty = 2.0)
  + `-t GONNET`: Gonnet160 substitution matrix (using path score, float type, gap open penalty = 22.0, gap extension penalty = 1.0)
- `-m`: Set the memory limit ratio. Default is 0.8, which sets the memory limit to 0.8 * available RAM. Set the ratio to negative to disable the memory-bound strategy.
  + Memory-bound strategy is disabled on macOS due to some technical issues.

For example, you can compute the exact MSA on `twi_009` without memory-bound strategy using the following command:

```
./RAM-MSA -f "../data/twi_009" -m -1
```

Or compute the exact MSA on `twi_009` with Gonnet160 substitution matrix and compute the SP accuracy and TC accuracy with its reference alignment:

```
./RAM-MSA -f "../data/twi_009" -rf "../data/twi_009_ref.fasta" -t GONNET
```

## Other details

1. Output file: `anytime_results.txt`
2. The running time limit is set to 10 million seconds (about 4 months). The search terminates if the timie limit is exceeded.
3. The memory usage is tracked by `/proc/meminfo` and `/proc/self/statm` so only available on linux. Disabled memory-bounded extension by `#if defined(__APPLE__) && defined(__MACH__)`.

## Trouble shooting

- CMake Error: CMake 3.15 or higher is required.
  + Solution: In `CMakeLists.txt` line 9, change `VERSION 3.14` to `VERSION 3.13` or lower. Tested with 3.15 and 3.14.
