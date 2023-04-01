#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void Usage(char* prog_name);
double **generate_matrix(int n);

int main(int argc, char* argv[]){
    int  n, num_threads;

    if (argc != 2) Usage(argv[0]);
    n = strtoll(argv[1], NULL, 10);

    double **A = generate_matrix(n);
    
    printf("Matriz generada aleatoriamente:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.3f  ", A[i][j]);
        }
        printf("\n");
    }

    num_threads = omp_get_max_threads();
    printf("Using %d threads\n", num_threads);
                      
    double eigenvalues[2];
    double eigenvectors[2][2];
    #pragma omp parallel sections shared(A, eigenvalues, eigenvectors)
    { 
        #pragma omp section
        {
			double discriminant = sqrt(pow(A[0][0] + A[1][1], 2) - 4 * (A[0][0] * A[1][1] - A[0][1] * A[1][0]));
        }
    }
    return 0;
}

//la matriz debe ser cuadrada para eigenvalores e inversa por lo que solo se pide un valor para generar la matriz
void Usage(char* prog_name) {
   fprintf(stderr, "Usage: %s <n>\n", prog_name);
   fprintf(stderr, "   n: number of rows and columns in matrix\n");
}


double **generate_matrix(int n) {
    int i,j;
    double **A = (double **) malloc(n * sizeof(double *));
    #pragma omp parallel for private(i,j) shared(A)
    for ( i = 0; i < n; i++) {
        A[i] = (double *) malloc(n * sizeof(double));
        #pragma omp parallel for
        for ( j = 0; j < n; j++) {
            A[i][j] = ((double) rand() / (double) RAND_MAX) * 1000;
        }
    }
    return A;
}




