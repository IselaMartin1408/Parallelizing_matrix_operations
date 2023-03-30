#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(){
    int size = 3;
    double matrix[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};
    double inverse[size*size];

    printf("Original matrix:\n");
    print_matrix(matrix, size);

    inverse_matrix(matrix, inverse, size);

    printf("Inverse matrix:\n");
    print_matrix(inverse, size);

    return 0;
}


void print_matrix(double* matrix, int size){
    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            printf("%f ", matrix[i*size+j]);
        }
        printf("\n");
    }
}

void inverse_matrix(double* matrix, double* inverse, int size){
    #pragma omp parallel for
    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            inverse[i*size+j] = (i==j)? 1.0 : 0.0;
        }
    }

    for (int k=0; k<size; k++){
        double pivot = matrix[k*size+k];
        #pragma omp parallel for
        for (int j=0; j<size; j++){
            matrix[k*size+j] /= pivot;
            inverse[k*size+j] /= pivot;
        }
        #pragma omp parallel for
        for (int i=k+1; i<size; i++){
            double factor = matrix[i*size+k];
            for (int j=0; j<size; j++){
                matrix[i*size+j] -= factor * matrix[k*size+j];
                inverse[i*size+j] -= factor * inverse[k*size+j];
            }
        }
    }

    for (int k=size-1; k>=0; k--){
        #pragma omp parallel for
        for (int i=k-1; i>=0; i--){
            double factor = matrix[i*size+k];
            for (int j=0; j<size; j++){
                matrix[i*size+j] -= factor * matrix[k*size+j];
                inverse[i*size+j] -= factor * inverse[k*size+j];
            }
        }
    }
}

