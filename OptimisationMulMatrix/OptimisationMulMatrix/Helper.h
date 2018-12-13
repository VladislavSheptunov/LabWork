#ifndef HELPER_H
#define HELPER_H

//#include "MatrixAlgorithms.h"

#include "stdio.h"
#include "time.h"
#include "stdlib.h"
#include "math.h"
#include "conio.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#define MSG(msg)            printf("%s\n", ##msg)

#define EPS                 0.001
#define IS_NOT_EQUAL(x, y)  ( fabs( (x) - (y)) > (EPS) ? true : false )
#define RANDOM(range)       ( (double)rand() / RAND_MAX * ( range ) ) // (double)(rand())/RAND_MAX*(max - min) + min;
//#define SIGN(x)             ( (x) >= 0 ? 1 : -1)

#define SIZE_MATRIX(size)           ( (size) * (size) )
#define SIZE_MATRIX_IN_BYTES(size)  ( SIZE_MATRIX(size) * sizeof(double) )

#define CLEAR_MATRIX(ptr)             assert(ptr); free((ptr))

#define MATRIX_MEMCPY(dst, src, size) memcpy(dst, src, SIZE_MATRIX_IN_BYTES(size));     

static inline void swap(double *left, double *right) {
  double tmp;

  tmp = *left;
  *left = *right;
  *right = tmp;
}

static inline void swap_str(double *left, double *right, int str_size) {
  size_t size_byte = str_size * sizeof(double);

  double *tmp = (double*)calloc(str_size, sizeof(double));
  assert(tmp);

  memcpy(tmp, left, size_byte);
  memcpy(left, right, size_byte);
  memcpy(right, tmp, size_byte);

  free(tmp);
}

static inline double* init_matrix(int matrix_size) {
  double *matrix = (double*)calloc(SIZE_MATRIX(matrix_size), sizeof(double));
  assert(matrix);
  return matrix;
}

static inline void fill_matrix(double *matrix, int matrix_size, double range) {
  assert(matrix);
  srand((unsigned int)time(NULL));
  for (int i = 0; i < matrix_size; i++)
    for (int j = 0; j < matrix_size; j++)
      *(matrix + i * matrix_size + j) = RANDOM(range);
}

static inline void fill_unit_matrix(double *matrix, int matrix_size) {
  assert(matrix);
  for (int i = 0; i < matrix_size; i++) {
    for (int j = 0; j < matrix_size; j++) {
      *(matrix + i * matrix_size + j) = 0.0;
      if (i == j)
        *(matrix + i * matrix_size + j) = 1.0;
    }
  }
}

static inline void show_matrix(const double *matrix, int matrix_size, const char* msg) {
  assert(matrix);
  MSG(msg);
  for (int i = 0; i < matrix_size; i++) {
    for (int j = 0; j < matrix_size; j++)
      printf("%.3f  ", *(matrix + i * matrix_size + j));
    printf(" \n");
  }
  printf("\n");
}

static inline void print_result_experiment_to_file(
  const std::vector<int> &size_matrix, 
  const std::vector<double> &run_time_before_opt,
  const std::vector<double> &run_time_after_opt,
  const std::string &file_name
) {
  std::ofstream tbl_speed_up;

  //system("md \"../ExperimentalResults\"");

  tbl_speed_up.open("../ExperimentalResults/" + file_name + ".csv");

  tbl_speed_up << "Run time\\Matrix size;";
  for (auto it : size_matrix)
    tbl_speed_up << std::to_string(it) << ";";
  tbl_speed_up << "\n";

  tbl_speed_up << "Not optimization algorithm;";
  for (auto it : run_time_before_opt)
    tbl_speed_up << std::to_string(it) << ";";
  tbl_speed_up << "\n";

  tbl_speed_up << "Optimized algorithm;";
  for (auto it : run_time_after_opt)
    tbl_speed_up << std::to_string(it) << ";";
  tbl_speed_up << "\n";

  tbl_speed_up << "SpeedUp;";
  for (int it = 0; it < run_time_before_opt.size(); it++)
    tbl_speed_up << std::to_string(run_time_before_opt[it] / run_time_after_opt[it]) << ";";
  tbl_speed_up << "\n";

  tbl_speed_up.close();
}

extern void inverse_matrix(double *matrix, int matrix_size);
extern void multiplication_matrix(
  double *res_matrix,
  const double *left_matrix,
  const double *right_matrix,
  int matrix_size
);

static inline void check_result_mul_matrix(
  const double *res_matrix,
  const double *left_matrix,
  const double *right_matrix,
  int matrix_size, 
  int matrix_size_for_show
) {
  assert(res_matrix && left_matrix && right_matrix);

  const double *A = left_matrix;
  const double *B = right_matrix;
  double *inverseB = init_matrix(matrix_size);
  double *E = init_matrix(matrix_size);
  double *A_mul_E = inverseB;

  MATRIX_MEMCPY(inverseB, B, matrix_size);

  inverse_matrix(inverseB, matrix_size);

  multiplication_matrix(E, B, (const double *)inverseB, matrix_size);
  multiplication_matrix(A_mul_E, A, E, matrix_size);

  bool fcheck = true;
  for (int i = 0; i < matrix_size; i++) {
    for (int j = 0; j < matrix_size; j++) {
      if (IS_NOT_EQUAL(*(A + i * matrix_size + j), *(A_mul_E + i * matrix_size + j))) {
        fcheck = false;
        break;
      }
    }
  }

  if (matrix_size <= matrix_size_for_show) {
    show_matrix(A, matrix_size, "The left Matrix:");
    show_matrix(B, matrix_size, "The right Matrix:");

    if (!fcheck)
      show_matrix(res_matrix, matrix_size, "A * B:");
  }

  if (fcheck)
    MSG("SUCCESSFUL");
  else
    MSG("ERROR");

  CLEAR_MATRIX(E);
  CLEAR_MATRIX(inverseB);
}

#endif //HELPER_H
