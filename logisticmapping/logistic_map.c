
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv)
{
  // taking values from argv
  double rmin = atof(argv[1]);
  double rmax = atof(argv[2]);
  double rstep = atof(argv[3]);
  int Niter = atoi(argv[4]);
  int n;
  
  // declaring array x of 1000 elements
  double x[Niter];
  
  printf("r n x");
  
  // to use above array for N elements
  int rem = Niter%1000;
  int t = (Niter-rem)/1000;
  
  //to skip first 100 iterations
  int s=0;
  
  //creating file to store data for gnuplot
  FILE *fptr;
  fptr = fopen("myfile.txt", "w"); 
  
  //to iterate for each r
  for (double r=rmin ; r<=rmax ; r+=rstep)
  {
    n=0;
    x[0] = 0.5;   // first element is anything
    //n++;
    //printf("%f %d %f \n", r, n, x[0]);
    //fprintf(fptr, "%f %f \n", r, x[0]);

    for (int i=1; i<Niter ; i++)
    {
      if (s==0)
      {
        if (i==100)
          s=1;
        continue;
      }
      x[i]=r*x[i-1]*(1-x[i-1]);
      n++;
      printf("%f %d %f \n", r, n, x[i]);
      fprintf(fptr, "%f %f \n", r, x[i]);
    } // closing inner loop
  } // closing outer loop
  
  // closing file
  fclose(fptr);
  
  // piping to gnuplot
  FILE *pipe_gp = popen("gnuplot -p", "w");
  fputs("set terminal png \n ", pipe_gp);
  fputs("set output 'myfile.png' \n", pipe_gp);
  fputs("set xlabel 'x' \n", pipe_gp);
  fputs("set xrange [0:4] \n", pipe_gp);
  fputs("set yrange [0:1] \n", pipe_gp);
  fputs("plot 'myfile.txt' with points", pipe_gp);
  pclose(pipe_gp);
  return (0);
}
