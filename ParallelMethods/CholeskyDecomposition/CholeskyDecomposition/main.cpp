#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#define PARALLEL 1
#define SOFT_GRADER 0

#if !SOFT_GRADER
	#include <conio.h>
	#include <vector>
	#include <fstream>
	#define EPS 0.001
	#define MSG(msg) printf("%s\n", ##msg)
	#define IS_NOT_EQUAL(x, y) ( fabs( (x) - (y)) > (EPS) ? true : false )
	#define RANDOM(range) ( (double)rand() / RAND_MAX * ( range ) ) // (double)(rand())/RAND_MAX*(max - min) + min;
	#define SIGN(x) ( (x) >= 0 ? 1 : -1)
#endif

#define SIZE_MATRIX(size) ( (size) * (size) )
#define SIZE_MATRIX_IN_BYTES(size) ( SIZE_MATRIX(size) * sizeof(double) )
#define CLEAR(ptr) \
if ((ptr) != NULL) free((ptr))

static inline double* init_matrix(int matrix_size) {
	double *matrix = (double*)calloc(SIZE_MATRIX(matrix_size), sizeof(double));
	assert(matrix != NULL);
	return matrix;
}

static inline double* init_block(int block_size) {
	double *block = (double*)calloc(block_size, sizeof(double));
	assert(block != NULL);
	return block;
}

static inline void swap(double *left, double *right) {
	double tmp;

	tmp = *left;
	*left = *right;
	*right = tmp;
}

#if !SOFT_GRADER
static inline void matrix_multiplication(double *res_matrix, const double *left_matrix, const double *right_matrix, int matrix_size) {
#pragma omp parallel for shared(left_matrix, right_matrix, res_matrix) private(i, j, k) if (PARALLEL)
	for (int i = 0; i < matrix_size; i++) {
		for (int j = 0; j < matrix_size; j++) {
			*(res_matrix + i * matrix_size + j) = 0;
			for (int k = 0; k < matrix_size; k++) {
				*(res_matrix + i * matrix_size + j) += (*(left_matrix + i * matrix_size + k) *  *(right_matrix + k * matrix_size + j));
			}
		}
	}
}

static inline void matrix_transposition(double *matrix, int matrix_size) {
#pragma omp parallel for shared(matrix, right_matrix) if (PARALLEL)
	for (int i = 0; i < matrix_size; i++) {
		for (int j = i; j < matrix_size; j++) {
			if (i != j)
				swap((matrix + i * matrix_size + j), (matrix + j * matrix_size + i));
		}
	}
}

static inline void fill_square_matrix(double *matrix, int matrix_size, double range) {
	assert(matrix != NULL);
	srand((unsigned int)time(NULL));
	for (int i = 0; i < matrix_size; i++)
		for (int j = 0; j < matrix_size; j++)
			*(matrix + i * matrix_size + j) = RANDOM(range);
}

static inline void gen_square_symmetric_positive_definite_matrix(double *matrix, int matrix_size, double range) {
	double *M = init_matrix(matrix_size);
	double *M_t = init_matrix(matrix_size);

	fill_square_matrix(M_t, matrix_size, range);
	memcpy(M, M_t, SIZE_MATRIX_IN_BYTES(matrix_size));
	matrix_transposition(M_t, matrix_size);
	matrix_multiplication(matrix, M, M_t, matrix_size);

	CLEAR(M);
	CLEAR(M_t);
}

static inline void show_matrix(const double *matrix, int matrix_size, const char* msg) {
	assert(matrix != NULL);
	MSG(msg);
	for (int i = 0; i < matrix_size; i++) {
		for (int j = 0; j < matrix_size; j++)
			printf("%.3f  ", *(matrix + i * matrix_size + j));
		printf(" \n");
	}
	printf("\n");
}

static inline void print_to_file(const std::vector<int> &size_matrix, const std::vector<int> &thread, const std::vector<double> &speedup) {
	system("del ./TableSpeedUp.csv");

	std::ofstream TableSpeedUp;
	//TableSpeedUp.open("TableSpeedUp.csv");
	TableSpeedUp.open("TableEfficiency.csv");
	TableSpeedUp << "Count thread\\Matrix size;";

	for (auto it : size_matrix)
		TableSpeedUp << it << ";";

	TableSpeedUp << "\n";

	for (int it = 0; it < thread.size(); it++) {
		TableSpeedUp << thread[it] << ";";
		for (int i = 0; i < size_matrix.size(); i++)
			TableSpeedUp << speedup[it + i*thread.size()] <<";";
		TableSpeedUp << "\n";
	}

	TableSpeedUp.close();
}

static inline void check_result_cholesky_decomposition(double *A, double *L, int  matrix_size, int matrix_size_for_show) {
	double *L_t = init_matrix(matrix_size);
	double *L_mul_L_t = init_matrix(matrix_size);

	memcpy(L_t, L, SIZE_MATRIX_IN_BYTES(matrix_size));
	matrix_transposition(L_t, matrix_size);
	matrix_multiplication(L_mul_L_t, L, L_t, matrix_size);

	CLEAR(L_t);

	bool fcheck = true;
	for (int i = 0; i < matrix_size; i++) {
		for (int j = 0; j < matrix_size; j++) {
			if (IS_NOT_EQUAL(*(A + i * matrix_size + j), *(L_mul_L_t + i * matrix_size + j))) {
				fcheck = false;
				break;
			}
		}
	}

	if (matrix_size <= matrix_size_for_show) {
		show_matrix(A, matrix_size, "The source Matrix:");
		show_matrix(L, matrix_size, "The lower triangular matrix:");
		// Result
		if (!fcheck)
			show_matrix(L_mul_L_t, matrix_size, "L * L^t:");
	}

	if (fcheck)
		MSG("SUCCESSFUL");
	else
		MSG("ERROR");

	CLEAR(L_mul_L_t);
}
#endif

/// ===================================== HELPER ===================================== ///

static inline void matrix_multiplication_block_left(double *res_matrix,
	const double *left_matrix,
	const double *right_matrix,
	int left_size_1, int left_size_2,
	int right_size_1, int right_size_2,
	int begin_i_left, int begin_j_left, int total_left_size)
{
	assert(left_size_2 == right_size_1);
#pragma omp parallel for shared(res_matrix, left_matrix, right_matrix) if (PARALLEL)
	for (int i = 0; i < left_size_1; i++) {
		for (int j = 0; j < right_size_2; j++) {
			*(res_matrix + i * right_size_2 + j) = 0;
			for (int k = 0; k < left_size_2; k++)
				*(res_matrix + i * right_size_2 + j) +=
				(*(left_matrix + (begin_i_left + i)*total_left_size + begin_j_left + k) * *(right_matrix + k * right_size_2 + j));
		}
	}
}

static inline void matrix_multiplication_block_res(double *res_matrix,
	const double *left_matrix,
	const double *right_matrix,
	int left_size_1, int left_size_2,
	int right_size_1, int right_size_2,
	int begin_i_left, int begin_j_left, int total_left_size,
	int begin_i_res, int begin_j_res, int total_res_size)
{
	assert(left_size_2 == right_size_1);
#pragma omp parallel for shared(res_matrix, left_matrix, right_matrix) if (PARALLEL)
	for (int i = 0; i < left_size_1; i++) {
		for (int j = 0; j < right_size_2; j++) {
			*(res_matrix + (begin_i_res + i)*total_res_size + begin_j_res + j) = 0;
			for (int k = 0; k < left_size_2; k++)
				*(res_matrix + (begin_i_res + i)*total_res_size + begin_j_res + j) +=
				(*(left_matrix + (begin_i_left + i)*total_left_size + begin_j_left + k) * *(right_matrix + k * right_size_2 + j));
		}
	}
}

static inline void matrix_subtraction_block(double *res_matrix,
	const double *left_matrix,
	const double *right_matrix,
	int size_1, int size_2,
	int begin_i, int begin_j, int total_size)
{
#pragma omp parallel for shared(res_matrix, left_matrix, right_matrix) if (PARALLEL)
	for (int i = 0; i < size_1; i++) {
		for (int j = 0; j < size_2; j++)
			*(res_matrix + i * size_1 + j) = *(left_matrix + (begin_i + i)*total_size + begin_j + j) - *(right_matrix + i * size_1 + j);
	}
}

static inline void matrix_transposition_block(double *matrix, const double *block,
	int size_1, int size_2,
	int begin_i, int begin_j, int total_size)
{
#pragma omp parallel for shared(matrix, block) if (PARALLEL)
	for (int i = 0; i < size_1; i++) {
		for (int j = 0; j < size_2; j++)
			*(matrix + j * size_1 + i) = *(block + (begin_i + i)*total_size + begin_j + j);
	}
}

static inline void get_cofactor(const double *matrix, double *tmp, int p, int q, int matrix_size) {
	int i = 0, j = 0;
	for (int row = 0; row < matrix_size; row++) {
		for (int col = 0; col < matrix_size; col++) {
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

static double determinant_matrix(double *matrix, int matrix_size) {
	if (matrix_size == 1)
		return *matrix;

	double *tmp = init_matrix(matrix_size);
	double det = 0;
	int sign = 1;

	for (int k = 0; k < matrix_size; k++) {
		get_cofactor(matrix, tmp, 0, k, matrix_size);
		det += sign * *(matrix + k) * determinant_matrix(tmp, matrix_size - 1);
		sign = -sign;
	}

	CLEAR(tmp);
	return det;
}

static inline void adjoint_matrix(double *A, double *adj, int matrix_size) {
	if (matrix_size == 1) {
		*adj = 1;
		return;
	}

#pragma omp parallel for shared(A, adj) if (PARALLEL)
	for (int i = 0; i < matrix_size; i++) {
		int sign = 1;
		double *tmp = init_matrix(matrix_size);
		for (int j = 0; j < matrix_size; j++) {
			get_cofactor(A, tmp, i, j, matrix_size);
			sign = ((i + j) % 2 == 0) ? 1 : -1;
			*(adj + j * matrix_size + i) = (sign)*(determinant_matrix(tmp, matrix_size - 1));
		}
		CLEAR(tmp);
	}

}

static inline void inverse_matrix(double *A, double *inverse, int matrix_size) {
	double det = determinant_matrix(A, matrix_size);
	assert(det != 0);
	double *adj = init_matrix(matrix_size);
	adjoint_matrix(A, adj, matrix_size);
#pragma omp parallel for shared(inverse, adj) if (PARALLEL)
	for (int i = 0; i < matrix_size; i++)
		for (int j = 0; j < matrix_size; j++)
			*(inverse + i * matrix_size + j) = *(adj + i * matrix_size + j) / det;
	CLEAR(adj);
}

/// ================================================================================== ///

/// ===================================== ALGORITHM ===================================== ///

#define BLOCK_SIZE  8

static double   *L11T_preallocated = NULL,
*L11T_inverse_preallocated = NULL;

static inline void Cholesky_Decomposition_line(const double *A, double *L, int n, int offsetof_size_A,
	int offsetof_col_L, int offsetof_row_L, int offsetof_size_L)
{
	/// https://ru.wikipedia.org/wiki/Разложение_Холецкого
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			sum = 0.0;
#pragma omp parallel for shared(L, A) reduction(+:sum) if (PARALLEL)
			for (int k = 0; k < j; k++)
				sum += (*(L + (i + offsetof_col_L) * offsetof_size_L + (k + offsetof_row_L)) *
					*(L + (j + offsetof_col_L) * offsetof_size_L + (k + offsetof_row_L)));
			*(L + (j + offsetof_row_L) + (i + offsetof_col_L) * offsetof_size_L) =
				((*(A + i * offsetof_size_A + j) - sum) / *(L + (j + offsetof_row_L) + (j + offsetof_col_L) * offsetof_size_L));
		}
		sum = 0.0;
#pragma omp parallel for shared(L, A) reduction(+:sum) if (PARALLEL)
		for (int k = 0; k < i; k++)
			sum += pow(*(L + (k + offsetof_row_L) + (i + offsetof_col_L) * offsetof_size_L), 2);
		*(L + (i + offsetof_row_L) + (i + offsetof_col_L) * offsetof_size_L) = sqrt(*(A + i * offsetof_size_A + i) - sum);
	}
}

static inline void Cholesky_Decomposition_Second_Iteration(double *A21, double *L11, double *L21, int n, int r,
	int begin_i_L11, int begin_j_L11, int total_len_L11,
	int begin_i_L21, int begin_j_L21, int total_len_L21,
	int begin_i_A21, int begin_j_A21, int total_len_A21)
{
	double	*L11T = L11T_preallocated,
		*L11T_inverse = L11T_inverse_preallocated;

	matrix_transposition_block(L11T, L11, r, r, begin_i_L11, begin_j_L11, total_len_L11);

	inverse_matrix(L11T, L11T_inverse, r);
	matrix_multiplication_block_res(L21, A21, L11T_inverse, (n - r), r, r, r,
		begin_i_A21, begin_j_A21, total_len_A21,
		begin_i_L21, begin_j_L21, total_len_L21);
}

static inline void Cholesky_Decomposition_Third_Iteration(double *A22_red, double *A22, double *L21, int n, int r,
	int begin_i_L21, int begin_j_L21, int total_len_L21,
	int begin_i_A22, int begin_j_A22, int total_len_A22)
{
	double	*L21_L21T = NULL,
		*L21T = init_block(r * (n - r));
#pragma omp parallel sections
	{
#pragma omp section
		{
			L21_L21T = init_matrix(n - r);
		}
#pragma omp section
		{
			matrix_transposition_block(L21T, L21, (n - r), r, begin_i_L21, begin_j_L21, total_len_L21);
		}
	}

	matrix_multiplication_block_left(L21_L21T,
		L21,
		L21T,
		(n - r), r,
		r, (n - r),
		begin_i_L21, begin_j_L21, total_len_L21
	);

#pragma omp parallel sections
	{
#pragma omp section
		{
			CLEAR(L21T);
		}
#pragma omp section
		{
			matrix_subtraction_block(A22_red, A22, L21_L21T, (n - r), (n - r), begin_i_A22, begin_j_A22, total_len_A22);
		}
	}
	CLEAR(L21_L21T);
}

static inline void Cholesky_Decomposition_Recursive(double *A, double *L, int n,
	double *A_full, int n_full,
	int L_begin_i, int L_begin_j)
{
	if (n <= BLOCK_SIZE) {
		Cholesky_Decomposition_line(A, L, n, n, L_begin_i, L_begin_j, n_full);
		CLEAR(A);
		return;
	}

	Cholesky_Decomposition_line(A, L, BLOCK_SIZE, n, L_begin_i, L_begin_j, n_full);

	Cholesky_Decomposition_Second_Iteration(A, L, L, n, BLOCK_SIZE,
		L_begin_i, L_begin_j, n_full,
		L_begin_i + BLOCK_SIZE, L_begin_j, n_full,
		BLOCK_SIZE, 0, n
	);


	double *A22_red = init_matrix(n - BLOCK_SIZE);

	Cholesky_Decomposition_Third_Iteration(A22_red, A, L, n, BLOCK_SIZE,
		L_begin_i + BLOCK_SIZE, L_begin_j + 0, n_full,
		BLOCK_SIZE, BLOCK_SIZE, n
	);

	CLEAR(A);

	Cholesky_Decomposition_Recursive(A22_red, L, n - BLOCK_SIZE, A_full, n_full,
		L_begin_i + BLOCK_SIZE, L_begin_j + BLOCK_SIZE);
}

void Cholesky_Decomposition(double *A, double *L, int n) {
	if (n <= BLOCK_SIZE) {
		Cholesky_Decomposition_line(A, L, n, n, 0, 0, n);
		return;
	}
#pragma omp parallel sections
	{
#pragma omp section
		{
			L11T_preallocated = init_matrix(BLOCK_SIZE);
			L11T_inverse_preallocated = init_matrix(BLOCK_SIZE);
		}
#pragma omp section
		{
			Cholesky_Decomposition_line(A, L, BLOCK_SIZE, n, 0, 0, n);
		}
	}

	Cholesky_Decomposition_Second_Iteration(A, L, L, n, BLOCK_SIZE,
		0, 0, n,
		BLOCK_SIZE, 0, n,
		BLOCK_SIZE, 0, n
	);

	double *A22_red = init_matrix(n - BLOCK_SIZE);

	Cholesky_Decomposition_Third_Iteration(A22_red, A, L, n, BLOCK_SIZE,
		BLOCK_SIZE, 0, n,
		BLOCK_SIZE, BLOCK_SIZE, n
	);

	Cholesky_Decomposition_Recursive(A22_red, L, n - BLOCK_SIZE, A, n, BLOCK_SIZE, BLOCK_SIZE);

	CLEAR(L11T_preallocated);
	CLEAR(L11T_inverse_preallocated);
}

/// ===================================================================================== ///

#if !SOFT_GRADER
int main(int argc, char** argv) {
	// Starting from the debug
	std::vector<int> arr_matrix_size = {250, 500, 1000, 1500, 2000, 2500, 3000};
	//std::vector<int> arr_matrix_size = { 5, 10 };
	std::vector<int> count_thread;
	std::vector<double> speedUp;
	int max_threads = omp_get_max_threads();
	for (int i = 1; i < max_threads + 1; i++) 
		count_thread.push_back(i);

	double *start_time = init_block(max_threads + 1);
	double *diff_time = init_block(max_threads + 1);
	double *speedup = init_block(max_threads + 1);
	double *efficiency = init_block(max_threads + 1);
	//--

	for (auto it : arr_matrix_size) {
		printf("\n========== MATRIX SIZE: %dx%d =========\n\n", it, it);
		omp_set_num_threads(1);
		double *matrix = init_matrix(it);
		gen_square_symmetric_positive_definite_matrix(matrix, it, 3.0);
		double *L = init_matrix(it);

		start_time[0] = omp_get_wtime();
		Cholesky_Decomposition_line(matrix, L, it, it, 0, 0, it);
		diff_time[0] = omp_get_wtime() - start_time[0];
		speedup[0] = 1;
		efficiency[0] = speedup[0] / 1;

		printf("The result of the sequential algorithm:\n");
		printf("Run time:   %.5f \n", diff_time[0]);
		printf("Speedup:    %.5f \n", speedup[0]);
		printf("Efficiency: %.5f \n\n", efficiency[0]);

		check_result_cholesky_decomposition(matrix, L, it, 10);

		for (int i = 1; i < max_threads + 1; i++) {
			omp_set_num_threads(i);
			memset(L, 0, SIZE_MATRIX_IN_BYTES(it));
			start_time[i] = omp_get_wtime();
			Cholesky_Decomposition(matrix, L, it);
			diff_time[i] = omp_get_wtime() - start_time[i];
			speedup[i] = diff_time[0] / diff_time[i];
			efficiency[i] = speedup[i] / i;
			printf("\nThe result of the block parallel algorithm (number of threads %d):\n", i);
			printf("Run time:   %.5f \n", diff_time[i]);
			printf("Speedup:    %.5f \n", speedup[i]);
			printf("Efficiency: %.5f \n\n\n", efficiency[i]);
			//speedUp.push_back(speedup[i]);
			speedUp.push_back(efficiency[i]);
			check_result_cholesky_decomposition(matrix, L, it, 10);
		}

		CLEAR(matrix);
		CLEAR(L);
	}

	print_to_file(arr_matrix_size, count_thread, speedUp);

	CLEAR(start_time);
	CLEAR(diff_time);
	CLEAR(speedup);
	CLEAR(efficiency);


	MSG("To end, press any key");
	_getch();

	return 0;
}
#endif