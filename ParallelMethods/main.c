#include "helper.h"
///*---------------------- Cholesky Decomposition -----------------------*///
#include "cholesky_decomposition.h"
///*--------------------------------------------------------------------*///

int main(int argc, char** argv){
    // Starting from the debug
    size_t size = 20;
    double *matrix = NULL;
    double *L = NULL;

    //double testMatrix[] = { 4.0, 12.0, -16.0, 12.0, 37.0, -43.0, -16.0, -43.0, 98.0 };

//    // Starting from the console
//    if (argc < 2){
//        printf("Invalid number of arguments!\n");
//        exit(-1);
//    }
//
//    if (RANGE(1,atoi(argv[1]),100000)){
//        size = atoi(argv[1]);
//    }else{
//        printf("Invalid matrix size!\n");
//        exit(-2);
//    }

    matrix = init_matrix(size);
    L = init_matrix(size);
    gen_square_symmetric_positive_definite_matrix(matrix, size, 5.0);
    //Cholesky_Decomposition_line(matrix, L, (int)size, (int)size, 0, 0, (int)size);
    Cholesky_Decomposition(matrix, L, (int)size);
    check_result_cholesky_decomposition(matrix, L, size, 10);

    CLEAR(matrix);
    CLEAR(L);

    printf ("Press Enter to continue...");
    getchar ();

    return 0;
}
//
//void Cholesky_Decomposition_line(const double *A, double *L, int n, int offsetof_size_A,
//                                 int offsetof_col_L, int offsetof_row_L, int offsetof_size_L)
//{
//    /// https://ru.wikipedia.org/wiki/Разложение_Холецкого
//    //memset(L, 0, SIZE_MATRIX(n));
//    double sum = 0;
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < i; j++) {
//            sum = 0;
//#pragma omp parallel for reduction(+:sum)
//            for (int k = 0; k < j; k++)
//                sum += ( *(L + (i + offsetof_col_L) * offsetof_size_L + (k + offsetof_row_L)) *
//                        *(L + (j + offsetof_col_L) * offsetof_size_L + (k + offsetof_row_L)) );
//            *(L + (j + offsetof_row_L) + (i + offsetof_col_L) * offsetof_size_L) =
//                    (( *(A + i*offsetof_size_A + j) - sum ) / *(L + (j + offsetof_row_L) + (j + offsetof_col_L) * offsetof_size_L ));
//        }
//        sum = 0;
//#pragma omp parallel for reduction(+:sum)
//        for (int k = 0; k < i; k++)
//            sum += pow(*(L + (k + offsetof_row_L) + (i + offsetof_col_L) * offsetof_size_L), 2);
//        *(L + (i + offsetof_row_L) + (i + offsetof_col_L) * offsetof_size_L) = sqrt(*(A + i*offsetof_size_A + i) - sum);
//    }
//}
//
//static inline void Cholesky_Decomposition_Second_Iteration(double *A21, double *L11, double *L21, int n, int r,
//                                                           int begin_i_L11, int begin_j_L11, int total_len_L11,
//                                                           int begin_i_L21, int begin_j_L21, int total_len_L21,
//                                                           int begin_i_A21, int begin_j_A21, int total_len_A21)
//{
//    double *L11T = L11T_preallocated, *L11T_inverse = L11T_inverse_preallocated;
//
//    matrix_transposition_block(L11T, L11, (size_t)r, (size_t)r, (size_t)begin_i_L11, (size_t)begin_j_L11, (size_t)total_len_L11);
//
//    inverse_matrix(L11T, L11T_inverse, (size_t)r);
//    matrix_multiplication_block_res(L21, A21, L11T_inverse, (n - r), r, r, r,
//                                    begin_i_A21, begin_j_A21, total_len_A21,
//                                    begin_i_L21, begin_j_L21, total_len_L21);
//}
//
//static void Cholesky_Decomposition_Third_Iteration(double *A22_red, double *A22, double *L21, int n, int r,
//                                                   int begin_i_L21, int begin_j_L21, int total_len_L21,
//                                                   int begin_i_A22, int begin_j_A22, int total_len_A22)
//{
//    double *L21_L21T = NULL,
//            *L21T = init_matrix((size_t)(r * (n - r)));
//#pragma omp parallel sections
//    {
//#pragma omp section
//        {
//            L21_L21T = init_matrix((size_t)((n - r) * (n - r)));
//        }
//#pragma omp section
//        {
//            matrix_transposition_block(L21T, L21,
//                                       (size_t)(n - r), (size_t)r,
//                                       (size_t)begin_i_L21, (size_t)begin_j_L21, (size_t)total_len_L21
//            );
//        }
//    }
//
//    matrix_multiplication_block_left(L21_L21T,
//                                     L21,
//                                     L21T,
//                                     (n - r), r,
//                                     r, (n - r),
//                                     begin_i_L21, begin_j_L21, total_len_L21
//    );
//
//#pragma omp parallel sections
//    {
//#pragma omp section
//        {
//            CLEAR(L21T);
//        }
//#pragma omp section
//        {
//            matrix_subtraction_block(A22_red,
//                                     A22,
//                                     L21_L21T,
//                                     (size_t)(n - r), (size_t)(n - r),
//                                     begin_i_A22, begin_j_A22, total_len_A22
//            );
//        }
//    }
//    CLEAR(L21_L21T);
//}
//
//static void Cholesky_Decomposition_Recursive(double *A, double *L, int n,
//                                             double *A_full, int n_full,
//                                             int L_begin_i, int L_begin_j)
//{
//    if (n <= BLOCK_SIZE) {
//        Cholesky_Decomposition_line(A, L, n, n, L_begin_i, L_begin_j, n_full);
//        CLEAR(A);
//        return;
//    }
//
//    Cholesky_Decomposition_line(A, L, BLOCK_SIZE, n, L_begin_i, L_begin_j, n_full);
//
//    Cholesky_Decomposition_Second_Iteration(A, L, L, n, BLOCK_SIZE,
//                                          L_begin_i, L_begin_j, n_full,
//                                          L_begin_i + BLOCK_SIZE, L_begin_j, n_full,
//                                          BLOCK_SIZE, 0, n
//    );
//
//
//    double *A22_red = init_matrix((size_t)(n - BLOCK_SIZE) * (n - BLOCK_SIZE));
//
//    Cholesky_Decomposition_Third_Iteration(A22_red, A, L, n, BLOCK_SIZE,
//                                           L_begin_i + BLOCK_SIZE, L_begin_j + 0, n_full,
//                                           BLOCK_SIZE, BLOCK_SIZE, n
//    );
//
//    CLEAR(A);
//
//    Cholesky_Decomposition_Recursive(A22_red, L, n - BLOCK_SIZE, A_full, n_full,
//                                     L_begin_i + BLOCK_SIZE, L_begin_j + BLOCK_SIZE);
//}
//
//void Cholesky_Decomposition(double *A, double *L, int n){
//    if (n <= BLOCK_SIZE) {
//        Cholesky_Decomposition_line(A, L, n, n, 0, 0, n);
//        return;
//    }
//#pragma omp parallel sections
//    {
//#pragma omp section
//        {
//            L11T_preallocated = init_matrix(BLOCK_SIZE);
//            L11T_inverse_preallocated = init_matrix(BLOCK_SIZE);
//        }
//#pragma omp section
//        {
//            Cholesky_Decomposition_line(A, L, BLOCK_SIZE, n, 0, 0, n);
//        }
//    }
//
//    Cholesky_Decomposition_Second_Iteration(A, L, L, n, BLOCK_SIZE,
//                                          0, 0, n,
//                                          BLOCK_SIZE, 0, n,
//                                          BLOCK_SIZE, 0, n
//    );
//
//    double *A22_red = init_matrix((size_t)(n - BLOCK_SIZE));
//
//    Cholesky_Decomposition_Third_Iteration(A22_red, A, L, n, BLOCK_SIZE,
//                                       BLOCK_SIZE, 0, n,
//                                       BLOCK_SIZE, BLOCK_SIZE, n
//    );
//
//
//    Cholesky_Decomposition_Recursive(A22_red, L, n - BLOCK_SIZE, A, n, BLOCK_SIZE, BLOCK_SIZE);
//
//    CLEAR(L11T_preallocated);
//    CLEAR(L11T_inverse_preallocated);
//}