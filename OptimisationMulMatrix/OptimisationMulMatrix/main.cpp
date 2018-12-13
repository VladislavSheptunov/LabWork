#include "windows.h"

#include "MatrixAlgorithms.h"

#define OPTIMAL 0

#define SIZE              1000
#define SIZE_MATRIX_SHOW  5
#define RANDOM_RANGE      5.0

int main(int argc, char** argv) {
  double start = 0;
  double runtime;
  double *left_matrix = init_matrix(SIZE);
  double *rigth_matrix = init_matrix(SIZE);
  double *res_matrix = init_matrix(SIZE);

  fill_matrix(left_matrix, SIZE, RANDOM_RANGE);
  fill_matrix(rigth_matrix, SIZE, RANDOM_RANGE);

  start = omp_get_wtime(); 
  {
    multiplication_matrix(res_matrix, left_matrix, rigth_matrix, SIZE);
  }
  runtime = omp_get_wtime() - start;  printf("Matrix multiplication algorithm. Not optimization:   %.5f seconds\n", runtime);
  check_result_mul_matrix(res_matrix, left_matrix, rigth_matrix, SIZE, SIZE_MATRIX_SHOW);

  start = omp_get_wtime();
  {
    optimal_multiplication_matrix(res_matrix, left_matrix, rigth_matrix, SIZE);
  }
  runtime = omp_get_wtime() - start;  printf("Matrix multiplication algorithm. Optimization:   %.5f seconds\n", runtime);
  check_result_mul_matrix(res_matrix, left_matrix, rigth_matrix, SIZE, SIZE_MATRIX_SHOW);

  CLEAR_MATRIX(left_matrix);
  CLEAR_MATRIX(rigth_matrix);
  CLEAR_MATRIX(res_matrix);

  MSG("To end, press any key");
  _getch();

  return EXIT_SUCCESS;
}