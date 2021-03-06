#include <vector>
#include "time.h"
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <omp.h>

#define PARALLEL 1
#define SOFT_GRADER 0

#if !SOFT_GRADER
	#include "conio.h"

	#define MSG(msg) printf("%s\n", ##msg)
	#define IS_EQUAL(x, y, eps) ( fabs( (x) - (y)) < (eps) ? true : false )
	#define RANDOM(range) ( (double)(rand()) / RAND_MAX * ( 2 * (range) )  - (range) )
#endif

#define SIZE_VECTOR(size) ( (size) * sizeof(double) )
#define CLEAR(ptr) \
if ((ptr) != NULL) free((ptr))

#if !SOFT_GRADER
// Why store the size of column (b), if the matrix is square (n x n)?
typedef struct _CRSMatrix {
	int n;						// Number of rows in the matrix
	int m;						// Number of columns in the matrix 
	int nz;						// Number of non-zero elements in a sparse matrix 
	std::vector<double> val;	// Array of matrix values by rows 
	std::vector<int> colIndex;	// Array of column numbers
	std::vector<int> rowPtr;	// Array of row start indices
}CRSMatrix;
#endif

static inline void init_crs_matrix(CRSMatrix &crs_matrix, int size, int notNull) {
	crs_matrix.n = crs_matrix.m = size;
	crs_matrix.nz = notNull;
	crs_matrix.val.resize(notNull);
	crs_matrix.colIndex.resize(notNull);
	crs_matrix.rowPtr.resize(size + 1);
}

static inline void clear_crs_matrix(CRSMatrix &crs_matrix) {
	crs_matrix.n = crs_matrix.m = crs_matrix.nz = 0;
	crs_matrix.val.clear();
	std::vector<double>().swap(crs_matrix.val);
	crs_matrix.colIndex.clear();
	std::vector<int>().swap(crs_matrix.colIndex);
	crs_matrix.rowPtr.clear();
	std::vector<int>().swap(crs_matrix.rowPtr);
}

static inline double* init_vector(int vector_size) {
	double *vector = (double*)malloc(SIZE_VECTOR(vector_size));
	assert(vector != nullptr);
	return vector;
}

#if !SOFT_GRADER
static inline void set_vector(double *vector, int size, double range) {
	srand((unsigned int)time(NULL));
	for (int i = 0; i < size; i++)
		vector[i] = RANDOM(range);
}

static inline void generate_crs_matrix(CRSMatrix &crs_matrix, int size, int count_nz_in_row, double range) {
	int f;
	init_crs_matrix(crs_matrix, size, count_nz_in_row * size);
	srand((unsigned int)time(NULL));
	for (int i = 0; i < crs_matrix.n; i++) {
		for (int j = 0; j < count_nz_in_row; j++) {
			do {
				crs_matrix.colIndex[i * count_nz_in_row + j] = rand() % crs_matrix.n;
				f = 0;
				for (int k = 0; k < j; k++)
					if (crs_matrix.colIndex[i * count_nz_in_row + j] ==
						crs_matrix.colIndex[i * count_nz_in_row + k])
						f = 1;
			} while (f == 1);
		}

		for (int j = 0; j < count_nz_in_row - 1; j++)
			for (int k = 0; k < count_nz_in_row - 1; k++)
				if (crs_matrix.colIndex[i * count_nz_in_row + k] >
					crs_matrix.colIndex[i * count_nz_in_row + k + 1])
				{
					std::swap(crs_matrix.colIndex[i * count_nz_in_row + k],
						crs_matrix.colIndex[i * count_nz_in_row + k + 1]);
				}
	}

	int c;
	for (int i = 0; i < count_nz_in_row * crs_matrix.n; i++) {
		crs_matrix.val[i] = RANDOM(range);
		c = 0;
		for (int j = 0; j <= crs_matrix.n; j++)
		{
			crs_matrix.rowPtr[j] = c;
			c += count_nz_in_row;
		}
	}
}

static double get_elem(const CRSMatrix &crs_matrix, int i, int j) {
	assert((i < crs_matrix.n) && (j < crs_matrix.n));
	double find_elem = 0.0;
	int j1 = crs_matrix.rowPtr[i];
	int j2 = crs_matrix.rowPtr[i + 1];
	for (int k = j1; k < j2; ++k) {
		if (crs_matrix.colIndex[k] == j) {
			find_elem = crs_matrix.val[k];
			break;
		}
	}
	return find_elem;
}
#endif

static inline void transpose_CRSMatrix(const CRSMatrix &crs_matrix, CRSMatrix &crs_matrix_t) {
	std::fill(crs_matrix_t.rowPtr.begin(), crs_matrix_t.rowPtr.end(), 0);

	for (int i = 0; i < crs_matrix.nz; i++)
		crs_matrix_t.rowPtr[crs_matrix.colIndex[i] + 1]++;

	int S = 0, tmp = 0;
	/// @TODO You can parallelize
	for (int i = 1; i <= crs_matrix.n; i++) {
		tmp = crs_matrix_t.rowPtr[i];
		crs_matrix_t.rowPtr[i] = S;
		S += tmp;
	}

	int j1 = 0, j2 = 0, column = 0;
	for (int i = 0; i < crs_matrix_t.n; i++) {
		j1 = crs_matrix.rowPtr[i];
		j2 = crs_matrix.rowPtr[i + 1];
		column = i;
		for (int j = j1; j < j2; j++) {
			crs_matrix_t.val[crs_matrix_t.rowPtr[crs_matrix.colIndex[j] + 1]] = crs_matrix.val[j];
			crs_matrix_t.colIndex[crs_matrix_t.rowPtr[crs_matrix.colIndex[j] + 1]] = column;
			crs_matrix_t.rowPtr[crs_matrix.colIndex[j] + 1]++;
		}
	}
}

static inline void mul_CRSMatrix_on_vector(const CRSMatrix &crs_matrix, const double *vector, double *result_vector) {
	int j1 = 0, j2 = 0;
	memset(result_vector, 0, SIZE_VECTOR(crs_matrix.n));
#pragma omp parallel for firstprivate(j1, j2) if (PARALLEL)
	for (int i = 0; i < crs_matrix.n; ++i) {
		j1 = crs_matrix.rowPtr[i];
		j2 = crs_matrix.rowPtr[i + 1];
		for (int j = j1; j < j2; ++j)
			result_vector[i] += crs_matrix.val[j] * vector[crs_matrix.colIndex[j]];
	}
}

static inline double mul_vector_on_vector(double *vector_1, double *vector_2, int size_vector) {
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum) if(PARALLEL)
	for (int i = 0; i < size_vector; i++)
		sum += vector_1[i] * vector_2[i];
	return sum;
}

void SLE_Solver_CRS_BICG(CRSMatrix &A, double *b, double eps, int max_iter, double *x, int &count) {
	// -- To speed up the calculations
	CRSMatrix At; 
	init_crs_matrix(At, A.n, A.nz);
	transpose_CRSMatrix(A, At);
	// --

	// -- Arrays for storing discrepancies current and next approximation
	double *R = init_vector(A.n);
	double *biR = init_vector(A.n);
	double *nR = init_vector(A.n);
	double *nbiR = init_vector(A.n);
	// --

	// -- Arrays for storing the current and next vector method step directions
	double *P = init_vector(A.n);
	double *biP = init_vector(A.n);
	double *nP = init_vector(A.n);
	double *nbiP = init_vector(A.n);
	// --

	// -- Arrays for storing the product of a matrix by a vector directions and bi-conjugate to it
	double *multAP = init_vector(A.n);
	double *multAtbiP = init_vector(A.n);
	// --

	// -- Coefficients of computational formulas
	double alfa, beta;
	// --

	// --Numerator and denominator of beta and alfa coefficients
	double numerator, denominator;
	// --

	// -- Variables to calculate the accuracy of the current approximation
	double check, norm;
	// --
	norm = sqrt(mul_vector_on_vector(b, b, A.n));

	int i; // Counter

	// -- Initialization
	std::fill_n(x, A.n, 1.0); // The initial approximation
	mul_CRSMatrix_on_vector(A, x, multAP);
	
	for (i = 0; i < A.n; i++)
		R[i] = biR[i] = P[i] = biP[i] = b[i] - multAP[i];

	// --

	// ----------------- ALGORITHM ----------------- //
	for (count = 0; count < max_iter; count++) {

		mul_CRSMatrix_on_vector(A, P, multAP);
		mul_CRSMatrix_on_vector(At, biP, multAtbiP);
		numerator = mul_vector_on_vector(biR, R, A.n);
		denominator = mul_vector_on_vector(biP, multAP, A.n);
		alfa = numerator / denominator;

#pragma omp parallel for shared(nR, R, multAP, nbiR, biR, multAtbiP) if (PARALLEL)
		for (i = 0; i < A.n; i++) {
			nR[i] = R[i] - alfa * multAP[i];
			nbiR[i] = biR[i] - alfa * multAtbiP[i];
		}

		denominator = numerator;
		numerator = mul_vector_on_vector(nbiR, nR, A.n);
		beta = numerator / denominator;

#pragma omp parallel for shared(nP, nR, P, nbiP, nbiR, biP, beta) if (PARALLEL)
		for (i = 0; i < A.n; i++) {
			nP[i] = nR[i] + beta * P[i];
			nbiP[i] = nbiR[i] + beta * biP[i];
		}

		check = sqrt(mul_vector_on_vector(R, R, A.n)) / norm;
		if (check < eps)
			break;

#pragma omp parallel for private(i) shared(x, P, alfa) if (PARALLEL)
		for (i = 0; i < A.n; i++) {
			x[i] += alfa * P[i];
		}

		std::swap(R, nR);
		std::swap(P, nP);
		std::swap(biR, nbiR);
		std::swap(biP, nbiP);
	}
	// -------------------------------------------- //
	
	clear_crs_matrix(At);
	CLEAR(R);
	CLEAR(biR);
	CLEAR(nR);
	CLEAR(nbiR);
	CLEAR(P);
	CLEAR(biP);
	CLEAR(nP);
	CLEAR(nbiP);
	CLEAR(multAP);
	CLEAR(multAtbiP);
}

#if !SOFT_GRADER
static inline void mul_CRSMatrix_on_vector_line(const CRSMatrix &crs_matrix, const double *vector, double *result_vector) {
	int j1 = 0, j2 = 0;
	memset(result_vector, 0, SIZE_VECTOR(crs_matrix.n));
	for (int i = 0; i < crs_matrix.n; ++i) {
		j1 = crs_matrix.rowPtr[i];
		j2 = crs_matrix.rowPtr[i + 1];
		for (int j = j1; j < j2; ++j)
			result_vector[i] += crs_matrix.val[j] * vector[crs_matrix.colIndex[j]];
	}
}

static inline double mul_vector_on_vector_line(double *vector_1, double *vector_2, int size_vector) {
	double sum = 0.0;
	for (int i = 0; i < size_vector; i++)
		sum += vector_1[i] * vector_2[i];
	return sum;
}

void SLE_Solver_CRS_BICG_line(CRSMatrix &A, double *b, double eps, int max_iter, double *x, int &count) {
	// -- To speed up the calculations
	CRSMatrix At;
	init_crs_matrix(At, A.n, A.nz);
	transpose_CRSMatrix(A, At);
	// --

	// -- Arrays for storing discrepancies current and next approximation
	double *R = init_vector(A.n);
	double *biR = init_vector(A.n);
	double *nR = init_vector(A.n);
	double *nbiR = init_vector(A.n);
	// --

	// -- Arrays for storing the current and next vector method step directions
	double *P = init_vector(A.n);
	double *biP = init_vector(A.n);
	double *nP = init_vector(A.n);
	double *nbiP = init_vector(A.n);
	// --

	// -- Arrays for storing the product of a matrix by a vector directions and bi-conjugate to it
	double *multAP = init_vector(A.n);
	double *multAtbiP = init_vector(A.n);
	// --

	// -- Coefficients of computational formulas
	double alfa, beta;
	// --

	// --Numerator and denominator of beta and alfa coefficients
	double numerator, denominator;
	// --

	// -- Variables to calculate the accuracy of the current approximation
	double check, norm;
	// --
	norm = sqrt(mul_vector_on_vector_line(b, b, A.n));

	int i; // Counter

		   // -- Initialization
	std::fill_n(x, A.n, 1.0); // The initial approximation
	mul_CRSMatrix_on_vector_line(A, x, multAP);

	for (i = 0; i < A.n; i++)
		R[i] = biR[i] = P[i] = biP[i] = b[i] - multAP[i];

	// --

	// ----------------- ALGORITHM ----------------- //
	for (count = 0; count < max_iter; count++) {

		mul_CRSMatrix_on_vector_line(A, P, multAP);
		mul_CRSMatrix_on_vector_line(At, biP, multAtbiP);
		numerator = mul_vector_on_vector_line(biR, R, A.n);
		denominator = mul_vector_on_vector_line(biP, multAP, A.n);
		alfa = numerator / denominator;

		for (i = 0; i < A.n; i++) {
			nR[i] = R[i] - alfa * multAP[i];
			nbiR[i] = biR[i] - alfa * multAtbiP[i];
		}

		denominator = numerator;
		numerator = mul_vector_on_vector_line(nbiR, nR, A.n);
		beta = numerator / denominator;

		for (i = 0; i < A.n; i++) {
			nP[i] = nR[i] + beta * P[i];
			nbiP[i] = nbiR[i] + beta * biP[i];
		}

		check = sqrt(mul_vector_on_vector_line(R, R, A.n)) / norm;
		if (check < eps)
			break;

		for (i = 0; i < A.n; i++) {
			x[i] += alfa * P[i];
		}

		std::swap(R, nR);
		std::swap(P, nP);
		std::swap(biR, nbiR);
		std::swap(biP, nbiP);
	}
	// -------------------------------------------- //

	clear_crs_matrix(At);
	CLEAR(R);
	CLEAR(biR);
	CLEAR(nR);
	CLEAR(nbiR);
	CLEAR(P);
	CLEAR(biP);
	CLEAR(nP);
	CLEAR(nbiP);
	CLEAR(multAP);
	CLEAR(multAtbiP);
}

static inline void show_crs_matrix(const CRSMatrix &crs_matrix, const char* msg) {
	MSG(msg);
	for (int i = 0; i < crs_matrix.m; i++) {
		for (int j = 0; j < crs_matrix.n; j++)
			printf("%.3f  ", get_elem(crs_matrix, i, j));
		printf("\n");
	}
	printf("\n");
}

static inline void show_vector(double *vector, int size, const char* msg) {
	MSG(msg);
	for (int i = 0; i < size; i++)
		printf("%.3f\n", vector[i]);
	printf("\n");
}

static inline void check_result(const CRSMatrix &crs_matrix, double *b, double *x, double eps) {
	double *check_vector = init_vector(crs_matrix.n);
	mul_CRSMatrix_on_vector(crs_matrix, x, check_vector);

	bool fcheck = true;
	for (int i = 0; i < crs_matrix.n; i++) {
		if (!IS_EQUAL(check_vector[i], b[i], eps)) {
			fcheck = false;
			break;
		}
	}

	if (crs_matrix.n <= 8) {
		show_crs_matrix(crs_matrix, "The sourse Matrix:");
		show_vector(b, crs_matrix.n, "Vector b:");
		show_vector(x, crs_matrix.n, "Vector x:");
		if (!fcheck)
			show_vector(check_vector, crs_matrix.n, "Check vector (A*x):");
	}

	if (fcheck)
		MSG("SUCCESSFUL");
	else
		MSG("ERROR");

	CLEAR(check_vector);
}

int main(int argc, char **argv) {
	int size = 8,
		max_iter = 10,
		count = 10;
	double eps = 0.001;

	
	CRSMatrix A;
	/// @TODO
	// For the correct operation of the algorithm, the matrix A must be nondegenerate, that is, detA! = 0.
	// In this case, the genius algorithm of the matrix is random. A variant of generating a degenerate matrix is possible.
	// Need to refine this algorithm
	generate_crs_matrix(A, size, 2, 2.0);

	double *b = init_vector(A.n);
	double *x = init_vector(A.n);
	set_vector(b, A.n, 2.0);
	
	clock_t start = clock();
	double wall_timer = omp_get_wtime();

	SLE_Solver_CRS_BICG(A, b, eps, max_iter, x, count);

	printf("TIME CLOCK = %f\n", (clock() - start) / 1000.0);
	printf("TIME WTIME = %f\n", (omp_get_wtime() - wall_timer));

	check_result(A, b, x, eps);

	clear_crs_matrix(A);
	CLEAR(b);
	CLEAR(x);

	MSG("To end, press any key");
	_getch();
	
	return 0;
}
#endif