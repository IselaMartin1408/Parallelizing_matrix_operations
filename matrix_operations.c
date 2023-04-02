#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

const int RMAX = 100;

void Usage(char* prog_name);
void Get_args(int argc, char* argv[], int* n_p, char* g_i_p);
void Print_matrix(double **A, int n, char* title);
void Generate_matrix(double **A, int n);
void Read_matrix(double **A, int n);
//double **matrix_inverse(double **A, int n);

int main(int argc, char* argv[]){
    int  n, num_threads;
    char g_i;
    double **A;
    
   Get_args(argc, argv, &n, &g_i);
   A = (double **) malloc(n * sizeof(double *));
    if (g_i == 'g') { //Hace una matriz random
        Generate_matrix(A, n);
   } else { //Pide al usuario que ingrese los datos de la matriz
        Read_matrix(A, n);
   }
   Print_matrix(A,n, "Matriz");//ocupen esta funcion para imprimir su matriz 

/*  
    double **A_inv = matrix_inverse(A, n);
    
    printf("Matriz inversa:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.3f  ", A_inv[i][j]);
        }
        printf("\n");
    }

    double **A_trans = transpose_matrix(A, n);
    
    printf("Matriz transpuesta:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.3f  ", A_trans[i][j]);
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
    }*/
    return 0;
}

/*-----------------------------------------------------------------
 * Function:  Usage
 * Purpose:    Resumen de cómo ejecutar un programa
 */
void Usage(char* prog_name) {
   fprintf(stderr, "Usage: %s <n>\n", prog_name);
   fprintf(stderr, "   n: number of rows and columns in matrix\n");
   fprintf(stderr, "  'g':  generate matrix using a random number generator\n");
   fprintf(stderr, "  'i':  user input matrix\n");
} /* Usage */


/*-----------------------------------------------------------------
 * Function:  Get_args
 * Purpose:   Obtener y comprobar los argumentos de la línea de comandos
 * In args:   argc, argv
 * Out args:  n_p, g_i_p
 */
void Get_args(int argc, char* argv[], int* n_p, char* g_i_p) {
   if (argc != 3 ) {
      Usage(argv[0]);
      exit(0);
   }
   *n_p = atoi(argv[1]);
   *g_i_p = argv[2][0];

   if (*n_p <= 0 || (*g_i_p != 'g' && *g_i_p != 'i') ) {
      Usage(argv[0]);
      exit(0);
   }
}  /* Get_args */

/*-----------------------------------------------------------------
 * Function:  Generate_matrix
 * Purpose:   Utilizar un generador de números aleatorios para generar los elementos de la matriz
 * In args:   n, A
 * Out args:  A
 */
void Generate_matrix(double **A, int n) {
    int i,j;

    #pragma omp parallel for private(i,j) shared(A)
    for ( i = 0; i < n; i++) {
        A[i] = (double *) malloc(n * sizeof(double));
        #pragma omp parallel for
        for ( j = 0; j < n; j++) {
            A[i][j] = ((double) rand() / (double) RAND_MAX) * 1000;
        }
    }
} /*Generate_matrix*/


/*-----------------------------------------------------------------
 * Function:  Print_matriz
 * Purpose:   Imprime los elementos de la matriz
 * In args:   A,n, titulo
 */

void Print_matrix(double **A, int n, char* title)
{
        printf("%s:\n", title);
        for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.3f  ", A[i][j]);
        }
        printf("\n");
    }
} /*Print_matriz*/


/*-----------------------------------------------------------------
 * Function:  Read_matriz
 * Purpose:   Lee los elementos ingresados por el usuario y los guarda en la matriz A
 * In args:   n,A
 * Out args:  A
 */
void Read_matrix(double **A, int n)
{
    printf("Ingrese los elementos de la matriz:\n");
    for (int i = 0; i < n; i++) {
        A[i] = (double *) malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
    }
} /*Read_matriz*/


/*

double **matrix_inverse(double **A, int n) {
    double **A_inv = (double **) malloc(n * sizeof(double *));
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        A_inv[i] = (double *) malloc(n * sizeof(double));
    }

    double det = 0;
    double **adj = (double **) malloc(n * sizeof (double *));
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        adj[i] = (double *) malloc(n * sizeof(double));
    }

    // Determinante
    for (int i = 0; i < n; i++) {
        det = det + (A[0][i] * (A[1][(i+1)%n] * A[2][(i+2)%n] - A[1][(i+2)%n] * A[2][(i+1)%n]));
    }

    if (det == 0) {
        printf("La matriz no es invertible.\n");
        return A_inv;
    }

    // Adjunta
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            adj[i][j] = pow(-1, i+j) * (A[(j+1)%n][(i+1)%n] * A[(j+2)%n][(i+2)%n] - A[(j+1)%n][(i+2)%n] * A[(j+2)%n][(i+1)%n]);
        }
    }

    // Matriz inversa
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_inv[i][j] = adj[i][j] / det;
        }
    }

    return A_inv;
}

}*/