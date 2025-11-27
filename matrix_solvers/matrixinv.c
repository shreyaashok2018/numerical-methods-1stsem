/* 
  Project : Matrix Inversion
  File : matrixinv.c
  Description : Uses gaussian elimination or lu decomposition to solve a particular set of linear equations.
  
  Implemented by : Shreya Ashok
  Collaborators : Deepa, Sumaiya 
  Semester : First Semester
  
  Usage : 
    gcc matrixinv.c -o matrixinv -lm
    ./matrixinv <filename> <N> <method> 
    
  Notes :
    Input file (<filename>) should contain an N x (N+1) augmented matrix: coefficients first, RHS in last column.
    N - the number of equations to be solved (hence, the file must have a minimum of N rows and N+1 columns)
    method - LU or Gaussian
    
  Example input file (for N = 3):
    2 1 -1 8
    -3 -1 2 -11
    -2 1 2 -3

  Assisted debugging with Gemini AI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void gauss(double **a,double *b,int N){
    double f;
    double x[N];
    double sum;
    double swap;

    //pivoting & elimination
    for (int k=0 ; k<=N-2 ; k++){
      int max = k;

    //finding row with largest element at kth column
    for (int i=k+1; i<=N-1 ; i++)
      max = (fabs(a[i][k])>fabs(a[max][k]))? i : max;
		        
    //swapping rows
    if ((a[max][k]!=0)&&(max!=k)){
      for(int i=0; i<=N-1; i++){
        swap = a[k][i];
	a[k][i] = a[max][i];
	a[max][i] = swap;
      }
      swap = b[k];
      b[k] = b[max];
      b[max] = swap;
    }
		        
    // elimination
    for (int i=k+1; i<=N-1 ; i++){
      f = a[i][k]/a[k][k];
      for (int j=k+1; j<=N-1 ; j++)
	a[i][j] -= f*a[k][j];
	b[i] -= f*b[k];
      }// for loop for i
    }
    
    //substitution
    for (int i=N-1; i>=0 ; i--){
      sum = b[i];
      for (int j=i+1 ; j<=N-1 ; j++)
        sum -= a[i][j]*x[j];
      x[i] = sum/a[i][i];
    }
    for (int i=0; i<N ; i++)
      printf("%lf \n",x[i]);        
}

void ludec(double **a,double *b,int N){
    double f;
    double x[N];
    double sum;
    double swap;
    double l[N][N];
    double d[N];
    
    for (int i=0; i<N ; i++){
      for (int j=0; j<N; j++)
          l[i][j] = (i==j)? 1:0;
    }    
    
    //pivoting & elimination
    for (int k=0 ; k<=N-2 ; k++){
		        
    // elimination
    for (int i=k+1; i<=N-1 ; i++){
      f = a[i][k]/a[k][k];
      for (int j=k+1; j<=N-1 ; j++)
	a[i][j] -= f*a[k][j];
      l[i][k] = f;
        }// for loop for i
    }
        
    //forward substitution
    for (int i=0; i<=N-1 ; i++){
      sum = b[i];
      for (int j=i-1 ; j>=0 ; j--)
        sum -= l[i][j]*d[j];
      d[i] = sum/l[i][i];
    }
        
    //backward substitution
    for (int i=N-1; i>=0 ; i--){
      sum = d[i];
      for (int j=i+1 ; j<=N-1 ; j++)
	sum -= a[i][j]*x[j];
      x[i] = sum/a[i][i];
    }
    for (int i=0; i<N ; i++)
      printf("%lf \n",x[i]);        
}

int main(int argc, char** argv){
    // Check command-line arguments
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <filename> <N> <Gaussian/LU>\n", argv[0]);
        return 1;
    }

    // Open input file
    char *name = argv[1];
    FILE *fp = fopen(name, "r");
    if (!fp) {
        perror("Failed to open input file");
        return 1;
    }
    
    int N = atoi(argv[2]);
    if (N<=0) {
      printf("N has to be a positive integer");
      exit(1);
      }
    
    char *method = argv[3];
	
    double **a = calloc(N, sizeof(double *));
    for (int i=0; i<N; i++) {
        a[i] = calloc(N, sizeof(double));
    }
    double *b = calloc(N, sizeof(double));
    
    //read file
    char line_buffer[8192]; //edit length as reqd

    //Read N rows
    for (int i = 0; i < N; i++) {
        
      //Read one full line from the file into the buffer
      if (fgets(line_buffer, sizeof(line_buffer), fp) == NULL) {
      
      // Error: reached the end of the file (or a read error occurred)
        fprintf(stderr, "Error or unexpected end-of-file at row %d\n", i);
        break; 
      }

      //Parse the N numbers from the line_buffer
      char *ptr = line_buffer; // Pointer to the current position in the buffer
      char *endptr;            // Pointer to store the end of the parsed number

      //Read N columns from the current line
      for (int j = 0; j < N; j++) {
              
        // Convert text to a double.
        double value = strtod(ptr, &endptr);

        // Check for parsing errors
        if (ptr == endptr) {
          // No number was found at the current position.
          fprintf(stderr, "Parse error on row %d, col %d. Line: %s\n", i, j, line_buffer);
          a[i][j] = 0.0; // Set a default value for remaining
        } 
        else a[i][j] = value;
              
        //Move our main pointer (`ptr`) to the end of number
        ptr = endptr;
      }
      b[i] = strtod(ptr, &endptr);
    }
    fclose(fp); 
	
    // code to read filename, open file and store in a dd array (a) & sd array (b)
	
    // check which method to use (copy strcomp from ode code)
    if (strcmp(method, "Gaussian" ) == 0) 
      gauss(a,b,N);
    else if (strcmp(method, "LU" ) == 0)
      ludec(a,b,N);
    else printf("Error: Unknown command. Enter Gaussian/LU \n");

    for (int i=0; i<N ; i++)
    free(a[i]);
    free(a);
    free(b);
    
    printf("\nSolution computed using %s method.\n", method);
    return 0;
}
