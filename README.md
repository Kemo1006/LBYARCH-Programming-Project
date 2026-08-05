# LBYARCH-Programming-Project

Write the kernel in (1) C program and (2) an x86-64 assembly language. The kernel is to perform DAXPY (A*X + Y) function.

## Input
- Scalar variable n (integer) contains the length of the vector;  Scalar variable A is a double-precision float. Vectors X, Y and Z are **double-precision float**.

    - Required to use functional scalar SIMD registers

    - Required to use functional scalar SIMD floating-point instructions

## Process: 
$Z[i] = A ⋅ X[i] + Y[i]$

**Example**:
    
    A --> 2.0

    x -> 1.0, 2.0, 3.0

    y -> 11.0, 12.0, 13.0

    (answer) z--> 13.0, 16.0, 19.0

**Output**: store result in vector Z.  Display the result of 1st ten elements of vector Z for all versions of kernel (i.e., C and x86-64).

## Note:

1. Write a C main program to call the kernels of the C version and x86-64 assembly language.

2. Time the kernel portion only.  

3. For each kernel version, time the process for vector size n = {220, 224, and  230}.  If 230 is impossible, you may reduce it to the point your machine can support (i.e.,  228 or 229).

4. You must run at least 30 times for each version to get the average execution time. 

5. For the data, you may initialize each vector and scalar variable with the same or different random value. 

6. You will need to check the correctness of your output.  Thus, if the C version is your "sanity check answer key," then the output of the x86-64 version has to be checked with the C version and output correspondingly (i.e., the x86-64 kernel output is correct, etc.).

7. Output in GitHub (make sure that I can access your Github):

a. Github readme containing the following (debug and release mode; C and x86-64):

    i. comparative execution time and short analysis of the performance of the kernels

    ii. Take a screenshot of the program output with the correctness check (C).

    iii. Take a screenshot of the program output, including the correctness check (x86-64).

    iv. short videos (5-10mins) showing your source code, compilation, and execution of the C and x86-64 program

b. Visual Studio project folder containing complete files (source code: C, x86-64, and all other required files) for others to load and execute your program.

## How to Build
1. Install NASM: https://www.nasm.us/
2. After installing, add NASM to environment variables
3. Windows Search **x64 Native Tools Command for VS**
4. To run
```
    cd /d "project folder"

    nasm -f win64 daxpy.asm -o daxpy.obj
    
    cl /O2 main.c daxpy.obj /Fe:daxpy_project.exe

    daxpy_project.exe
```