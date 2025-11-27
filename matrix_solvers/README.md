# Matrix Inversion – Gaussian Elimination & LU Decomposition

This program solves a system of **N linear equations** using either:

- **Gaussian Elimination** (with partial pivoting)  
- **LU Decomposition**

Both methods are implemented **from scratch** in C, without any external libraries.

---

## How to Compile & Run

```bash
gcc matrixinv.c -o matrixinv -lm
./matrixinv input.txt N Gaussian      # OR
./matrixinv input.txt N LU
```

# Input
input file must contain an augmented matrix with N rows and N+1 columns
sample input file input.txt is present in the folder - this works for N upto 10.

# Output
Each line of output will have one solution vector

# Collaboration
According to our professor's instructions, the project was done in groups of 3. Collaborated with Deepa and Sumaiya.

