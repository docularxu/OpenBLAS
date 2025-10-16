
# OpenBLAS Benchmark Guidelines

This guide demonstrates the process of benchmarking the OpenBLAS library on a Spacemit K1/M1 RISC-V processor.
It covers performance evaluation from baseline scalar execution to vectorized code utilizing the ZVL128B and ZVL256B extensions with different LMUL configurations.

## 1. Download
In this guide we download the latest release v0.3.30 to build the OpenBLAS instead of the `develop` branch to ensure the best reliability.
```bash
wget https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v0.3.30.tar.gz
tar zxvf v0.3.30.tar.gz
cd OpenBLAS-0.3.30
```

## 2. Compilation
### Building the default `TARGET`
The following table shows the builtin `TARGET` for respective compiling options.

|      TARGET       | float LMUL(M*N) | double LMUL(M*N) |  Type   |
|-------------------|-----------------|------------------|---------|
| `RISCV64_GENERIC` |       N/A       |        N/A       | Scalar  |
| `RISCV64_ZVL128B` |     2 (8*8)     |      4 (8*4)     | ZVL128B |
| `RISCV64_ZVL256B` |     1 (16*8)    |      1 (8*8)     | ZVL256B |

For compiling options other than the default `TARGET` (e.g., ZVL128B with LMUL=1), please refer to the next chapter.
```
# Set TARGET, for example, RISCV64_ZVL256B
TARGET="RISCV64_ZVL256B"

# 1. Compile the OpenBLAS library
make TARGET=$TARGET -j $(nproc)

# 2. Compile the benchmark utilities
cd benchmark
make TARGET=$TARGET -j $(nproc)
```
This will build the OpenBLAS lib and all the benchmark programs.
### Building custom `TARGET`
To build a specific TARGET (e.g., RISCV64_ZVL128B) with a custom LMUL value (e.g., LMUL=1), the simplest approach is to first generate the custom kernel and then modify the corresponding KERNEL.<TARGET> configuration file.

OpenBLAS provides the Python script kernel/riscv64/generate_kernel.py to generate the code template (referred to as a "kernel") for GEMM and GEMV, allowing customization of LMUL, M and N.

#### Example 1: Building ZVL128B with `LMUL=1` for Single Precision (float)
```bash
cd kernel/riscv64
python3 ./generate_kernel.py op=gemm param_precision=float complex=false LMUL=1 M=8 N=8 reg_width_bits=128 cpu=zvl128b output=sgemm_kernel_8x8_zvl128b_m1.c
```
Then modify the `kernel/riscv64/KERNEL.RISCV64_ZVL128B`, locating the line for `SGEMMKERNEL` and replace the default entry 
```
SGEMMKERNEL    =  sgemm_kernel_$(SGEMM_UNROLL_M)x$(SGEMM_UNROLL_N)_zvl128b.c
```
with the newly generated kernel file
```
SGEMMKERNEL    =  sgemm_kernel_8x8_zvl128b_m1.c
```
#### Example 2: Building ZVL256B with `LMUL=2` for Double Precision (double)
```bash
cd kernel/riscv64
python3 ./generate_kernel.py op=gemm param_precision=double complex=false LMUL=2 M=8 N=8 reg_width_bits=256 cpu=zvl256b output=dgemm_kernel_8x8_zvl256b_m2.c
```
Then modify the `kernel/riscv64/KERNEL.RISCV64_ZVL256B`, locating the line for `DGEMMKERNEL` and replace the default entry
```
DGEMMKERNEL    =  dgemm_kernel_$(DGEMM_UNROLL_M)x$(DGEMM_UNROLL_N)_zvl256b.c
```
with the newly generated kernel file
```
DGEMMKERNEL    =  dgemm_kernel_8x8_zvl256b_m2.c
```
------
Then follow the compilation steps presented in the last chapter to build the lib and benchmark.

> NOTE: Modifying the builtin `TARGET` configuration is the fastest way to customize the kernel. For a more elegant solution, like adding a completely new, separate `TARGET`, kindly follow the steps in [this pull request](https://github.com/OpenMathLib/OpenBLAS/pull/4355/files)

## 3. Benchmark
OpenBLAS has already provided utilities to test the performance of various operators. You can just use it right away.

For example, to test the GEMM performance in `float` precision. 

```bash
# OPENBLAS_LOOPS specify the loop count
# ./sgemm.goto <min matrix size> <max matrix size> <exponent>
# 
# In our modified version of gemm.c, the benchmark starts at the 
# minimum size and multiplies the size by the exponent in each step
# until the maximum size is reached or exceeded

$ OPENBLAS_LOOPS=100 ./sgemm.goto 64 256 2 > OPENBLAS_SGEMM_RISCV64_GENERIC_10_256-4096.csv
From :  64  To : 256 Step=2 : Transa=N : Transb=N
          SIZE                   Flops             Time
 M=  64, N=  64, K=  64 :     1489.98 MFlops   0.035188 sec
 M= 128, N= 128, K= 128 :    10431.79 MFlops   0.040207 sec
 M= 256, N= 256, K= 256 :    10413.86 MFlops   0.322209 sec
```

We strongly recommends setting the `OPENBLAS_LOOPS` to a value greater than 1 for much more stable result (the default value is 1, which leads to unstable result since the calculation time is too short).

In the modified versions of `gemm.c` and `gemv.c`, the raw performance data in CSV format is directed to `stdout`, while the real-time execution results are displayed on `stderr`. Therefore, you can easily redirect the CSV output to a file (e.g., `OPENBLAS_SGEMM_RISCV64_GENERIC_64-256.csv`) for subsequent analysis and plotting.
 
