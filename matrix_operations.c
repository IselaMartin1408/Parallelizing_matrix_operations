#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main(){
 
 int num_threads = omp_get_max_threads();
    printf("Using %d threads\n", num_threads);

    double A[2][2] = {{1.0, 2.0},
                      {3.0, 4.0}};
    double eigenvalues[2];
    double eigenvectors[2][2];
    #pragma omp parallel sections shared(A, eigenvalues, eigenvectors)
    { 
        #pragma omp section
        {
			double discriminant = sqrt(pow(A[0][0] + A[1][1], 2) - 4 * (A[0][0] * A[1][1] - A[0][1] * A[1][0]));
        }

    return 0;
}

