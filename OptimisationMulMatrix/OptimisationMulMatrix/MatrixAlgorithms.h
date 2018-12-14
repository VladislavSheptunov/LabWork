#ifndef _MATRIX_ALGORITHMS_H_
#define _MATRIX_ALGORITHMS_H_

#include "assert.h"
#include "omp.h"
#include "Helper.h"

static inline void transposition_matrix(double *matrix, int matrix_size) {
  assert(matrix);
  for (int i = 0; i < matrix_size; i++) {
    for (int j = i; j < matrix_size; j++) {
      if (i != j)
        swap((matrix + i * matrix_size + j), (matrix + j * matrix_size + i));
    }
  }
}

static inline void inverse_matrix(double *matrix, int matrix_size) {
  assert(matrix);

  double multi = 0;
  double *E = init_matrix(matrix_size);
  fill_unit_matrix(E, matrix_size);

  for (int k = 0; k < matrix_size; ++k) {
    if (fabs(*(matrix + k * matrix_size + k)) < 1e-8) {
      bool changed = false;

      for (int i = k + 1; i < matrix_size; ++i) {
        if (fabs(*(matrix + i * matrix_size + k)) < 1e-8) {
          swap_str((matrix + k), (matrix + i), matrix_size);
          swap_str((E + k), (E + i), matrix_size);
          changed = true;
          break;
        }
      }

      if (!changed)
        MSG("Matrix cannot inverse!");
    }

    double div = *(matrix + k * matrix_size + k);
    for (int j = 0; j < matrix_size; ++j) {
      *(matrix + k * matrix_size + j) /= div;
      *(E + k * matrix_size + j) /= div;
    }

#pragma omp parallel for firstprivate(k, multi)
    for (int i = k + 1; i < matrix_size; ++i) {
      multi = *(matrix + i * matrix_size + k);
      for (int j = 0; j < matrix_size; ++j) {
        *(matrix + i * matrix_size + j) -= multi * *(matrix + k * matrix_size + j);
        *(E + i * matrix_size + j) -= multi * *(E + k * matrix_size + j);
      }
    }
  }

  multi = 0;
  for (int k = matrix_size - 1; k > 0; --k) {
#pragma omp parallel for firstprivate(k, multi)
    for (int i = k - 1; i >= 0; --i) {
      multi = *(matrix + i * matrix_size + k);
      for (int j = 0; j < matrix_size; ++j) {
        *(matrix + i * matrix_size + j) -= multi * *(matrix + k * matrix_size + j);
        *(E + i * matrix_size + j) -= multi * *(E + k * matrix_size + j);
      }
    }
  }

  MEMCPY_MATRIX(matrix, E, matrix_size);

  DESTROY_MATRIX(E);
}

static inline void multiplication_matrix(
  double *res_matrix, 
  const double *left_matrix, 
  const double *right_matrix, 
  int matrix_size
) {

  assert(res_matrix && left_matrix && right_matrix);

  for (int i = 0; i < matrix_size; i++) {
    for (int j = 0; j < matrix_size; j++) {
      *(res_matrix + i * matrix_size + j) = 0;
      for (int k = 0; k < matrix_size; k++) {
        *(res_matrix + i * matrix_size + j) += (*(left_matrix + i * matrix_size + k) *  *(right_matrix + k * matrix_size + j));
      }
    }
  }
}

/* =============================================================================== /
                                Optimization steps
   1. Using Intel Commpiller
   2. Cache memory optimization
   3. Vectorizing Using AVX Intrinsics
   4. Paralleling with using OpenMP
/  =============================================================================== */

static inline void optimal_multiplication_matrix(
  double *res_matrix,
  const double *left_matrix,
  const double *right_matrix,
  int matrix_size
) {

  assert(res_matrix && left_matrix && right_matrix);

#ifdef COUNT_THREAD
  omp_set_num_threads(COUNT_THREAD);
#endif // COUNT_THREAD

#pragma omp parallel for if (PARALLELIZATION)
  for (int i = 0; i < matrix_size; ++i) {
    for (int k = 0; k < matrix_size; ++k) {
#ifdef VECTORIZATION
      __m256d v_res_matrix = _mm256_setzero_pd();
      for (int j = 0; j < matrix_size; j += 4) {
        __m256d v_left_matrix = _mm256_broadcast_sd(left_matrix + i * matrix_size + k);
        __m256d v_right_matrix = _mm256_loadu_pd(right_matrix + k * matrix_size + j);
        __m256d v_tmp_matrix = _mm256_mul_pd(v_left_matrix, v_right_matrix);
        v_res_matrix = _mm256_add_pd(v_res_matrix, v_tmp_matrix);
        _mm256_store_pd(res_matrix + i * matrix_size + j, v_res_matrix);
      }
#else
      for (int j = 0; j < matrix_size; ++j)
        *(res_matrix + i * matrix_size + j) += (*(left_matrix + i * matrix_size + k) *  *(right_matrix + k * matrix_size + j));
#endif
      }
    }
  }

#endif //_MATRIX_ALGORITHMS_H_
