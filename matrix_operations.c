#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

const int RMAX = 100;

void Usage(char* prog_name);
void Get_args(int argc, char* argv[], int* n_p, char* g_i_p);
void Print_matrix(double **A, int n, char* title);
void Generate_matrix(double **A, int n);
void Read_matrix(double **A, int n);
void EigenProgram(double A[2][2],double  eigenvalues[],double eigenvectors[][2]);
void Matrix_inverse(double **A, int n);
double Matrix_determinant(double **A, int n);
void Transpose(double **adj, int n);
void Mxm(double **A,double **C, int n);

int main(int argc, char* argv[]){
    int  n;
    char g_i;
    double **A;
    double **C;
    
	// variables eigenvalores
	double B[2][2] = {{1.0, 2.0},
                      {3.0, 4.0}};
    double eigenvalues[2];
    double eigenvectors[2][2];
	// fin variables eigenvalores
	
   Get_args(argc, argv, &n, &g_i);
   A = (double **) malloc(n * sizeof(double *));
   C = (double **) malloc(n * sizeof(double *));
    if (g_i == 'g') { //Hace una matriz random
        Generate_matrix(A, n);
        Generate_matrix(C, n);
   } else { //Pide al usuario que ingrese los datos de la matriz
        Read_matrix(A, n);
        Read_matrix(C, n);
   }
   Print_matrix(A,n, "Matriz A");
   Print_matrix(C,n, "Matriz B");


    Matrix_inverse(A, n);
	EigenProgram(B,eigenvalues,eigenvectors);
    Transpose(A,n); 
    Mxm(A, C, n);
	
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
    int i, j;
    unsigned int seed;

    // Generar una semilla aleatoria para cada fila
    for (i = 0; i < n; i++) {
        seed = (unsigned int) time(NULL) + i;
        srand(seed);
        
        // Generar elementos aleatorios para cada columna de la fila
        A[i] = (double*) malloc(n * sizeof(double));
        for (j = 0; j < n; j++) {
            A[i][j] = ((double) rand() / (double) RAND_MAX) * 100;
        }
    }
}  /*Generate_matrix*/


/*-----------------------------------------------------------------
 * Function:  Print_matriz
 * Purpose:   Imprime los elementos de la matriz
 * In args:   A,n, titulo
 */

void Print_matrix(double **A, int n, char* title)
{
        printf("\n-----%s-----\n", title);
        for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("[ %.3f ]", A[i][j]);
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



/*-----------------------------------------------------------------
 * Function:  Matrix_inverse
 * Purpose:   Calcula la inversa de la matriz
 * In args:   n,A
 * Out args:  A_inv
 */
void Matrix_inverse(double **A, int n) {
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
    if (n == 1) {
        det = A[0][0];
    } else {
        for (int i = 0; i < n; i++) {
            double **submatrix = (double **) malloc((n-1) * sizeof(double *));
            for (int j = 0; j < n-1; j++) {
                submatrix[j] = (double *) malloc((n-1) * sizeof(double));
            }
            for (int j = 1; j < n; j++) {
                int k = 0;
                for (int l = 0; l < n; l++) {
                    if (l != i) {
                        submatrix[j-1][k] = A[j][l];
                        k++;
                    }
                }
            }
            double sub_det = Matrix_determinant(submatrix, n-1);
            det += A[0][i] * pow(-1, i) * sub_det;
        }
    }

    if (det == 0) {
        printf("La matriz no es invertible.\n");
        return;
    }
    printf("\nDeterminante = %f\n", det);

    // Adjunta
    #pragma omp parallel for 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double **submatrix = (double **) malloc((n-1) * sizeof(double *));
            for (int k = 0; k < n-1; k++) {
                submatrix[k] = (double *) malloc((n-1) * sizeof(double));
            }
            int p = 0;
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    int q = 0;
                    for (int l = 0; l < n; l++) {
                        if (l != j) {
                            submatrix[p][q] = A[k][l];
                            q++;
                                               }
                }
                p++;
            }
        }
        double sub_det = Matrix_determinant(submatrix, n-1);
        adj[j][i] = pow(-1, i+j) * sub_det;
    }
}

        Print_matrix(adj,n, "Matriz adjunta");

// Matriz inversa
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        A_inv[i][j] = adj[i][j] / det;
    }
}

// Mostrar resultados
Print_matrix(A_inv,n, "Matriz inversa");

}/*Matrix_inverse*/


/*-----------------------------------------------------------------
 * Function:  EigenProgram
 * Purpose:   Calcula el eigenvalor e eigenvector de una matriz
 * In args:   A, eigenvalues, eigenvectors
 * Out args:  A_inv
 */
void EigenProgram(double A[2][2],double  eigenvalues[],double eigenvectors[][2]){
#pragma omp parallel sections shared(A, eigenvalues, eigenvectors)
    {
       
        #pragma omp section
        {
            double discriminant = sqrt(pow(A[0][0] + A[1][1], 2) - 4 * (A[0][0] * A[1][1] - A[0][1] * A[1][0]));
            eigenvalues[0] = (A[0][0] + A[1][1] + discriminant) / 2.0;
            eigenvalues[1] = (A[0][0] + A[1][1] - discriminant) / 2.0;
        }


        #pragma omp section
        {
            double lambda1 = eigenvalues[0];
            double lambda2 = eigenvalues[1];
            double c = A[1][0];
            double d = A[1][1];
            double norm1 = sqrt(pow(lambda1 - d, 2) + pow(c, 2));
            double norm2 = sqrt(pow(lambda2 - d, 2) + pow(c, 2));
            eigenvectors[0][0] = (lambda1 - d) / norm1;
            eigenvectors[0][1] = c / norm1;
            eigenvectors[1][0] = (lambda2 - d) / norm2;
            eigenvectors[1][1] = c / norm2;
        }
    }


    printf("\nEigenvalues: %f, %f\n", eigenvalues[0], eigenvalues[1]);
    printf("\n----Eigenvectors:----\n");
    printf("[ %f, %f ]\n", eigenvectors[0][0], eigenvectors[0][1]);
    printf("[ %f, %f ]\n", eigenvectors[1][0], eigenvectors[1][1]);

} /*EigenProgram*/


/*-----------------------------------------------------------------
 * Function:  Transpose
 * Purpose:   Calcula la transpuesta de una matriz
 * In args:   A, n
 * Out args:  A_inv
 */
void Transpose(double **A, int n) {
    int i, j;
    double **trans;
    trans = (double **) malloc(n * sizeof(double *));
for (int i = 0; i < n; i++) {
    trans[i] = (double *) malloc(n * sizeof(double));
}
     #pragma omp parallel for private(i,j) shared(trans)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            trans[j][i] = A[i][j];
        }
    }
    Print_matrix(trans,n, "Transpuesta");
}/*Transpose*/


/*-----------------------------------------------------------------
 * Function:  Matrix_determinant
 * Purpose:   Calcula la determinante de una matriz
 * In args:   A, n
 * Out args:  det
 */
double Matrix_determinant(double **A, int n) {
    double det = 0;
    if (n == 1) {
        det = A[0][0];
    } else {
        for (int i = 0; i < n; i++) {
            double **submatrix = (double **) malloc((n-1) * sizeof(double *));
            for (int j = 0; j < n-1; j++) {
                submatrix[j] = (double *) malloc((n-1) * sizeof(double));
            }
            for (int j = 1; j < n; j++) {
                int k = 0;
                for (int l = 0; l < n; l++) {
                    if (l != i) {
                        submatrix[j-1][k] = A[j][l];
                        k++;
                    }
                }
            }
            double sub_det = Matrix_determinant(submatrix, n-1);
            det += A[0][i] * pow(-1, i) * sub_det;
        }
    }
    return det;
}/*Matrix_determinant*/


/*-----------------------------------------------------------------
 * Function:  Mxm
 * Purpose:   Calcula la multipliacion de dos matrices
 * In args:   A, C, n
 * Out args:  MT
 */
void Mxm(double **A, double **C, int n)
{
    int i, j,k;
   double **MT;
   MT = (double **) malloc(n * sizeof(double *));
   for (int i = 0; i < n; i++) {
      MT[i] = (double *) malloc(n * sizeof(double));
   }
    // Inicializar la matriz C a cero
    #pragma omp parallel for private(j,k)
    for ( i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            MT[i][j] = 0;
        }
    }

    // Multiplicar las matrices A y B
    #pragma omp parallel for private(j,k)
    for ( i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for ( k = 0; k < n; k++) {
                MT[i][j] += A[i][k] * C[k][j];
            }
        }
    }
    Print_matrix(MT,n, "Multiplicacion");
} /*Mxm*/