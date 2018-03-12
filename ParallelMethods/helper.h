//
// Created by Vladislav on 11.03.2018.
//
#ifndef PARALLELMETHODS_HELPER_H
#define PARALLELMETHODS_HELPER_H

#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "math.h"
#include "stdbool.h"

// index = i*columns+j !

#define SIZE_MATRIX(size) ( (size) * (size) * sizeof(double) )
#define RANDOM(range) ( (double)(rand()) / RAND_MAX * ( 2 * (range) )  - (range) )
#define IS_EQUAL(x, y, eps) ( fabs( (x) - (y)) < (eps) ? true : false )
#define IS_NOT_EQUAL(x, y, eps) ( !IS_EQUAL(x, y, eps) )

#define IS_NULL(ptr) \
if ((ptr) == NULL) { printf("Pointer is NULL. File: %s, line: %d", __FILE__, __LINE__); exit(-1); }
#define CLEAR(ptr) \
if ((ptr) != NULL) free((ptr))

static inline double* init_matrix(size_t matrix_size){
    double *matrix = NULL;
    matrix = (double*)malloc(SIZE_MATRIX(matrix_size));
    IS_NULL(matrix);
    return matrix;
}

static inline void fill_square_matrix(double *matrix, size_t matrix_size, double range){
    IS_NULL(matrix);
    srand ((unsigned int)time(NULL));
    for(size_t i = 0; i < matrix_size; i++)
        for(size_t j = 0; j < matrix_size; j++)
            *(matrix + i * matrix_size + j) = RANDOM(range);
}

static inline void swap(double *left, double *right){
    double tmp;

    tmp = *left;
    *left = *right;
    *right = tmp;
}

static void matrix_multiplication(double *res_matrix,
                                  const double *left_matrix,
                                  const double *right_matrix,
                                  size_t matrix_size)
{
    size_t i,j,k;
#pragma omp parallel num_threads(4)
#pragma omp parallel for shared(left_matrix, right_matrix, res_matrix) private(i, j, k) if (matrix_size >= 4)
    for (i = 0; i < matrix_size; i++) {
        for (j = 0; j < matrix_size; j++) {
            *(res_matrix + i * matrix_size + j) = 0;
            for (k = 0; k < matrix_size; k++) {
                *(res_matrix + i * matrix_size + j) += ( *(left_matrix + i * matrix_size + k) *  *(right_matrix + k * matrix_size + j));
            }
        }
    }
}

static void matrix_multiplication_block_left(double *res_matrix,
                                            const double *left_matrix,
                                            const double *right_matrix,
                                            int left_size_1, int left_size_2,
                                            int right_size_1, int right_size_2,
                                            int begin_i_left, int begin_j_left, int total_left_size)
{
    if (left_size_2 != right_size_1){
        printf("Blocks not equals");
        exit(-3);
    }
#pragma omp parallel for
    for (size_t i = 0; i < left_size_1; i++) {
        for (size_t j = 0; j < right_size_2; j++) {
            *(res_matrix + i*right_size_2 + j) = 0;
            for (int k = 0; k < left_size_2; k++)
                *(res_matrix + i*right_size_2 + j) +=
                        (*(left_matrix + (begin_i_left + i)*total_left_size + begin_j_left + k) * *(right_matrix + k*right_size_2 + j));
        }
    }
}

static void matrix_multiplication_block_res(double *res_matrix,
                                            const double *left_matrix,
                                            const double *right_matrix,
                                            int left_size_1, int left_size_2,
                                            int right_size_1, int right_size_2,
                                            int begin_i_left, int begin_j_left, int total_left_size,
                                            int begin_i_res, int begin_j_res, int total_res_size)
{
    if (left_size_2 != right_size_1){
        printf("Blocks not equals");
        exit(-3);
    }
#pragma omp parallel for
    for (size_t i = 0; i < left_size_1; i++) {
        for (size_t j = 0; j < right_size_2; j++) {
            *(res_matrix + (begin_i_res + i)*total_res_size + begin_j_res + j) = 0;
            for (size_t k = 0; k < left_size_2; k++)
                *(res_matrix + (begin_i_res + i)*total_res_size + begin_j_res + j) +=
                (*(left_matrix + (begin_i_left + i)*total_left_size + begin_j_left + k) * *(right_matrix + k*right_size_2 + j));
        }
    }
}

static void matrix_subtraction_block(double *res_matrix,
                                     const double *left_matrix,
                                     const double *right_matrix_,
                                     size_t size_1, size_t size_2,
                                     int begin_i, int begin_j, int total_size)
{
#pragma omp parallel for
    for (size_t i = 0; i < size_1; i++) {
        for (size_t j = 0; j < size_2; j++)
            *(res_matrix + i*size_1 + j) = *(left_matrix + (begin_i + i)*total_size + begin_j + j) - *(right_matrix_ + i*size_1 + j);
    }
}

static inline void get_cofactor(const double *matrix, double *tmp, size_t p, size_t q, size_t matrix_size){
    size_t i = 0, j = 0;
    for (size_t row = 0; row < matrix_size; row++) {
        for (size_t col = 0; col < matrix_size; col++) {
            if (row != p && col != q) {
                *(tmp + i * (matrix_size - 1) + (j++)) = *(matrix + row * matrix_size + col);
                if (j == matrix_size - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

static double determinant_matrix(double *matrix, size_t matrix_size){
    if (matrix_size == 1)
        return *matrix;

    double *tmp = init_matrix(matrix_size);
    double det = 0;
    int sign = 1;

    for (size_t k = 0; k < matrix_size; k++) {
        get_cofactor(matrix, tmp, 0, k, matrix_size);
        det += sign * *(matrix + k) * determinant_matrix(tmp, matrix_size - 1);
        sign = -sign;
    }

    CLEAR(tmp);
    return det;
}

static inline void adjoint_matrix(double *A, double *adj, size_t matrix_size){
    if (matrix_size == 1) {
        *adj = 1;
        return;
    }
#pragma omp parallel for
    for (size_t i = 0; i < matrix_size; i++) {
        int sign = 1;
        double *tmp = init_matrix(matrix_size);
        for (size_t j = 0; j < matrix_size; j++) {
            get_cofactor(A, tmp, i, j, matrix_size);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            *(adj + j*matrix_size + i) = (sign)*(determinant_matrix(tmp, matrix_size - 1));
        }
        CLEAR(tmp);
    }
}

static inline void inverse_matrix(double *A, double *inverse, size_t matrix_size){
    double det = determinant_matrix(A, matrix_size);
    if (det == 0) {
        printf("Singular matrix, can't find its inverse");
        exit(-2);
    }

    double *adj = init_matrix(matrix_size);
    adjoint_matrix(A, adj, matrix_size);
#pragma omp parallel for
    for (size_t i = 0; i < matrix_size; i++)
        for (size_t j = 0; j < matrix_size; j++)
            *(inverse + i*matrix_size + j) = *(adj + i*matrix_size + j) / det;
    CLEAR(adj);
}

static inline void matrix_transposition(double *matrix, size_t matrix_size){
#pragma omp parallel for
    for(size_t i = 0; i < matrix_size; i++) {
        for (size_t j = i; j < matrix_size; j++) {
            if (i != j)
                swap((matrix + i * matrix_size + j), (matrix + j * matrix_size + i));
        }
    }
}

static inline void matrix_transposition_block(double *matrix, const double *block,
                                              size_t size_1, size_t size_2,
                                              size_t begin_i, size_t begin_j, size_t total_size)
{
#pragma omp parallel for
    for (size_t i = 0; i < size_1; i++) {
        for (size_t j = 0; j < size_2; j++)
            *(matrix + j*size_1 + i) = *(block + (begin_i + i)*total_size + begin_j + j);
    }
}

static inline void show_matrix(const double *matrix, size_t matrix_size){
    IS_NULL(matrix);
    for (size_t i = 0; i < matrix_size; i++) {
        for (size_t j = 0; j < matrix_size; j++)
            printf("%.3f  ", *(matrix + i * matrix_size +j));
        printf(" \n");
    }
    printf("\n");
}

static void gen_square_symmetric_positive_definite_matrix(double *matrix, size_t matrix_size, double range){
    double *M =     init_matrix(matrix_size);
    double *M_t =   init_matrix(matrix_size);

    fill_square_matrix(M_t, matrix_size, range);
    memcpy(M, M_t, SIZE_MATRIX(matrix_size));
    matrix_transposition(M_t, matrix_size);
    matrix_multiplication(matrix, M, M_t, matrix_size);

    CLEAR(M);
    CLEAR(M_t);
}

static void check_result_cholesky_decomposition(double *A, double *L, size_t matrix_size, size_t matrix_size_for_show){
    double *L_t =       init_matrix(matrix_size);
    double *L_mul_L_t = init_matrix(matrix_size);

    memcpy(L_t, L, SIZE_MATRIX(matrix_size));
    matrix_transposition(L_t, matrix_size);
    matrix_multiplication(L_mul_L_t, L, L_t, matrix_size);

    CLEAR(L_t);

    _Bool fcheck = true;
    for (size_t i = 0; i < matrix_size; i++) {
        for (size_t j = 0; j < matrix_size; j++) {
            if (IS_NOT_EQUAL(*(A + i * matrix_size +j), *(L_mul_L_t + i * matrix_size +j), 0.00001)){
                fcheck = false;
                break;
            }
        }
    }

    // Show Matrix
    if (matrix_size <= matrix_size_for_show){
        printf("The source Matrix:\n");
        show_matrix(A,matrix_size);
        printf("The lower triangular matrix:\n");
        show_matrix(L,matrix_size);
        // Result
        printf("L * L^t:\n");
        show_matrix(L_mul_L_t, matrix_size);
    }

    if (fcheck)
        printf("SUCCESSFUL\n");
    else
        printf("ERROR\n");

    CLEAR(L_mul_L_t);
}
#endif //PARALLELMETHODS_HELPER_H
