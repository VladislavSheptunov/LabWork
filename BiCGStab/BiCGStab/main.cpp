#include <vector>
#include "time.h"
#include <iostream>
#include <algorithm> 

#define EPS 0.001
#define SIZE_MATRIX(size) ( (size) * (size) * sizeof(double) )
#define SIZE_VECTOR(size) ( (size) * sizeof(double) )
#define RANDOM(range) ( (double)(rand()) / RAND_MAX * ( 2 * (range) )  - (range) )
#define CLEAR(ptr) \
if ((ptr) != NULL) free((ptr))
#define ERROR_MSG(msg) \
printf("%s. File: %s, line: %d\n", ##msg, __FILE__, __LINE__)
#define MSG(msg) \
printf("%s\n", ##msg)
#define IS_EQUAL(x, y, eps) ( fabs( (x) - (y)) < (eps) ? true : false )
#define IS_NOT_EQUAL(x, y, eps) ( !IS_EQUAL(x, y, eps) )

typedef struct _CRSMatrix {
	int n;						// „исло строк в матрице 
	int m;						// „исло столбцов в матрице 
	int nz;						// „исло ненулевых элементов в разреженной матрице 
	std::vector<double> val;	// ћассив значений матрицы по строкам 
	std::vector<int> colIndex;	// ћассив номеров столбцов 
	std::vector<int> rowPtr;	// ћассив индексов начала строк 
}CRSMatrix;

template <typename T>
static inline void swap(T *left, T *right) {
	T tmp;

	tmp = *left;
	*left = *right;
	*right = tmp;
}

static inline double* init_matrix(int matrix_size) {
	double *matrix = nullptr;
	matrix = (double*)malloc(SIZE_MATRIX(matrix_size));
	if (matrix == nullptr)
		ERROR_MSG("Pointer is NULL");
	return matrix;
}

static inline double* init_vector(int vector_size) {
	double *vector = nullptr;
	vector = (double*)malloc(SIZE_VECTOR(vector_size));
	if (vector == nullptr)
		ERROR_MSG("Pointer is NULL");
	return vector;
}

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
					swap(&crs_matrix.colIndex[i * count_nz_in_row + k],
						&crs_matrix.colIndex[i * count_nz_in_row + k + 1]);
				}
	}

	int c;
	for (int i = 0; i < count_nz_in_row * crs_matrix.n; i++) {
		crs_matrix.val[i] = RANDOM(range);
		c = 0;
		for (int j = 0; j <= crs_matrix.n ; j++)
		{
			crs_matrix.rowPtr[j] = c;
			c += count_nz_in_row;
		}
	}
}

static double get_elem(const CRSMatrix &crs_matrix, int i, int j) {
	if (crs_matrix.n < i || crs_matrix.n < j)
		ERROR_MSG("Error in index");
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

static inline void show_matrix(const double *matrix, int matrix_size) {
	if (matrix == nullptr)
		ERROR_MSG("Pointer is NULL");
	for (int i = 0; i < matrix_size; i++) {
		for (int j = 0; j < matrix_size; j++)
			printf("%.3f  ", *(matrix + i * matrix_size + j));
		printf(" \n");
	}
	printf("\n");
}

static inline void show_crs_matrix(const CRSMatrix &crs_matrix, const char* msg) {
	MSG(msg);
	for (int i = 0; i < crs_matrix.m; i++) {
		for (int j = 0; j < crs_matrix.n; j++) {
			if (j != 0)
				printf(" ");
			printf("%.3f", get_elem(crs_matrix, i, j));
		}

		if (i < crs_matrix.m) {
			printf("\n");
		}
	}
	printf("\n");
}

static inline void show_vector(double *vector, int size, const char* msg) {
	MSG(msg);
	for (int i = 0; i < size; i++)
		printf("%.3f\n", vector[i]);
	printf("\n");
}

static inline void mtrxcpy(double *matrix, const CRSMatrix &crs_matrix) {
	if (matrix == nullptr)
		ERROR_MSG("Pointer is NULL");
	int i1 = 0, 
		i2 = 0;
	for (int i = 0; i < crs_matrix.n; i++) {
		i1 = crs_matrix.rowPtr[i];
		i2 = crs_matrix.rowPtr[i + 1];
		for (int j = i1; j < i2; j++)
			matrix[i * crs_matrix.n + crs_matrix.colIndex[j]] = crs_matrix.val[j];
	}
}

static int crsmatrixcmp(const CRSMatrix &crs_matrix_1, const CRSMatrix &crs_matrix_2) {
	if (crs_matrix_1.n != crs_matrix_2.n)
		return 1;

	int ret = 0;

	double *tmp_matrix = init_matrix(crs_matrix_1.n);
	memset(tmp_matrix, 0, SIZE_MATRIX(crs_matrix_1.n));
	mtrxcpy(tmp_matrix, crs_matrix_1);

	int i1 = 0,
		i2 = 0;
	for (int i = 0; i < crs_matrix_2.n; ++i) {
		i1 = crs_matrix_2.rowPtr[i];
		i2 = crs_matrix_2.rowPtr[i + 1];
		for (int j = i1; j < i2; j++)
			tmp_matrix[i * crs_matrix_2.n + crs_matrix_2.colIndex[j]] -= crs_matrix_2.val[j];
	}

	for (int i = 0; i < crs_matrix_2.n; i++)
		for (int j = 0; j < crs_matrix_2.n; j++)
			if (tmp_matrix[i * crs_matrix_2.n + j] > EPS)
				ret = 1;

	CLEAR(tmp_matrix);

	return ret;
}

static inline void transpose_CRSMatrix(const CRSMatrix &crs_matrix, CRSMatrix &crs_matrix_t) {
	std::fill(crs_matrix_t.rowPtr.begin(), crs_matrix_t.rowPtr.end(), 0);
	for (int i = 0; i < crs_matrix.nz; i++)
		crs_matrix_t.rowPtr[crs_matrix.colIndex[i] + 1]++;

	int S = 0, tmp = 0;
	for (int i = 1; i <= crs_matrix.n; i++) {
		tmp = crs_matrix_t.rowPtr[i];
		crs_matrix_t.rowPtr[i] = S;
		S += tmp;
	}

	int j1 = 0, j2 = 0, 
		col = 0, row = 0, 
		IIndex = 0, RIndex = 0;
	//double val;
	for (int i = 0; i < crs_matrix_t.n; i++) {
		j1 = crs_matrix.rowPtr[i]; 
		j2 = crs_matrix.rowPtr[i + 1];
		col = i; // —толбец в AT - строка в ј
		for (int j = j1; j < j2; j++) {
			//val = crs_matrix.val[j];
			RIndex = crs_matrix.colIndex[j];
			IIndex = crs_matrix_t.rowPtr[RIndex + 1];
			//crs_matrix_t.val[IIndex] = val;
			crs_matrix_t.val[IIndex] = crs_matrix.val[j];
			crs_matrix_t.colIndex[IIndex] = col;
			crs_matrix_t.rowPtr[RIndex + 1]++;
		}
	}
}

static inline void mul_CRSMatrix_on_vector(const CRSMatrix &crs_matrix, const double *vector, double *result_vector) {
	int s = 0, k = 0;
	memset(result_vector, 0, SIZE_VECTOR(crs_matrix.n));
	for (int i = 0; i < crs_matrix.n; ++i) {
		s = crs_matrix.rowPtr[i];
		k = crs_matrix.rowPtr[i + 1];
		for (int j = s; j < k; ++j)
			result_vector[i] += crs_matrix.val[j] * vector[crs_matrix.colIndex[j]];
	}
}

static inline double mul_vector_on_vector(double *vector_1, double *vector_2, int size_vector) {
	double sum = 0.0;
	for (int i = 0; i < size_vector; i++)
		sum += vector_1[i] * vector_2[i];
	return sum;
}

static inline void mul_vector_on_scalar(double *vector, int size_vector, double scalar) {
	for (int i = 0; i < size_vector; i++)
		vector[i] *= scalar;
}

static inline void check_result(const CRSMatrix &crs_matrix, double *b, double *x, double eps) {
	double *check_vector = init_vector(crs_matrix.n);
	mul_CRSMatrix_on_vector(crs_matrix, x, check_vector);

	bool fcheck = true;
	for (int i = 0; i < crs_matrix.n; i++) {
		if (IS_NOT_EQUAL(check_vector[i], b[i], eps)) {
			fcheck = false;
			break;
		}
	}

	if (crs_matrix.n <= 5) {
		show_crs_matrix(crs_matrix, "CRS MATRIX A:");
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

void SLE_Solver_CRS_BICG(CRSMatrix &A, double *b, double eps, int max_iter, double *x, int &count) {
	// ƒл€ ускорени€ вычислений вычислим
	// транспонированную матрицу ј
	CRSMatrix At;
	init_crs_matrix(At, A.n, A.nz);

	transpose_CRSMatrix(A, At);

	// массивы дл€ хранени€ нев€зки
	// текущего и следующего приближени€
	double *R = init_vector(A.n);
	double *biR = init_vector(A.n);
	double *nR = init_vector(A.n);
	double *nbiR = init_vector(A.n);

	// массивы дл€ хранени€ текущего и следующего вектора
	// направлени€ шага метода
	double *P = init_vector(A.n);
	double *biP = init_vector(A.n);
	double *nP = init_vector(A.n);
	double *nbiP = init_vector(A.n);

	// указатель, дл€ смены указателей на вектора текущего
	// и следующего шага метода
	double * tmp;

	// массивы дл€ хранени€ произведени€ матрицы на вектор
	//направлени€ и бисопр€женный к нему
	double *multAP = init_vector(A.n);
	double *multAtbiP = init_vector(A.n);

	// beta и alfa - коэффициенты расчетных формул
	double alfa, beta;
	// числитель и знаменатель коэффициентов beta и alfa
	double numerator, denominator;
	// переменные дл€ вычислени€
	// точности текущего приближени€
	double check, norm;
	norm = sqrt(mul_vector_on_vector(b, b, A.n));

	//задание начального приближени€
	int i;
	//int n = A.N;
	memset(x, 1, SIZE_VECTOR(A.n));
	//инициализаци€ метода
	mul_CRSMatrix_on_vector(A, x, multAP);

	for (i = 0; i < A.n; i++)
		R[i] = biR[i] = P[i] = biP[i] = b[i] - multAP[i];

	// реализаци€ метода
	for (count = 0; count < max_iter; count++) {
		mul_CRSMatrix_on_vector(A, P, multAP);
		mul_CRSMatrix_on_vector(At, biP, multAtbiP);
		numerator = mul_vector_on_vector(biR, R, A.n);
		denominator = mul_vector_on_vector(biP, multAP, A.n);
		alfa = numerator / denominator;

		for (i = 0; i < A.n; i++)
			nR[i] = R[i] - alfa * multAP[i];

		for (i = 0; i < A.n; i++)
			nbiR[i] = biR[i] - alfa * multAtbiP[i];

		denominator = numerator;
		numerator = mul_vector_on_vector(nbiR, nR, A.n);
		beta = numerator / denominator;

		for (i = 0; i < A.n; i++)
			nP[i] = nR[i] + beta * P[i];

		for (i = 0; i < A.n; i++)
			nbiP[i] = nbiR[i] + beta * biP[i];

		// контроль достежени€ необходимой точности
		check = sqrt(mul_vector_on_vector(R, R, A.n)) / norm;
		if (check < eps)
			break;

		for (i = 0; i < A.n; i++)
			x[i] += alfa * P[i];

		// мен€ем массивы текущего и следующего шага местами
		tmp = R; R = nR; nR = tmp;
		tmp = P; P = nP; nP = tmp;
		tmp = biR; biR = nbiR; nbiR = tmp;
		tmp = biP; biP = nbiP; nbiP = tmp;
	}

	// освобождение пам€ти
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

int main(char **argv, int argc) {
	int size = 3,
		max_iter = 10,
		count = 10;
	double eps = 0.001;

	CRSMatrix A;

	generate_crs_matrix(A, size, 2, 2.0);

	double *b = init_vector(A.n);
	set_vector(b, A.n, 2.0);
	double *x = init_vector(A.n);

	SLE_Solver_CRS_BICG(A, b, eps, max_iter, x, count);

	check_result(A, b, x, eps);

	CLEAR(b);
	CLEAR(x);

	system("pause");
	return 0;
}