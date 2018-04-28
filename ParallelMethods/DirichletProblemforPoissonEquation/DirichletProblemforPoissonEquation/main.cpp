//#include <iostream>
//#include <omp.h>
//#include <iomanip>
//#include <vector>
//#include <cmath>
//#include <cstring>
//#include <memory>
//#include <numeric>
//#include <algorithm>

//
//#define SOFT_GRADER 0
//#define PARALLEL 1
//
//
//#define M_PI	3.14159265358979323846  
//#define EPS		1e-4f  
//
//typedef std::vector<std::vector<double>> grid;
//
//#if !SOFT_GRADER

#define HEAD_TASK_1 1
#define HEAD_TASK_2 0

#if HEAD_TASK_1
class heat_task {
	heat_task();
public:
	double X;	// ширина пластины
	double Y;	// длина пластины
	int n;		// размер сетки по x
	int m;		// размер сетки по y
	explicit heat_task(double _x, double _y, int _n, int _m) : X(_x), Y(_y), n(_n), m(_m) {}
	double left_condition(double y)		{ return 0.0 + y * y;	} // функция, задающая граничное условие при x = 0
	double right_condition(double y)	{ return 4.0 + y * y;	} // функция, задающая граничное условие при x = X
	double bottom_condition(double x)	{ return 0.0 + x * x;	} // функция, задающая граничное условие при y = 0
	double top_condition(double x)		{ return 9.0 + x * x;	} // функция, задающая граничное условие при y = Y
	double f(double x, double y)		{ return 4.0;			} // функция, задающая внешнее воздействие
};
#endif // HEAD_TASK_1

#if HEAD_TASK_2
class heat_task {
	heat_task();
public:
	double X;	// ширина пластины
	double Y;	// длина пластины
	int n;		// размер сетки по x
	int m;		// размер сетки по y
	explicit heat_task(double _x, double _y, int _n, int _m) : X(_x), Y(_y), n(_n), m(_m) {}
	double left_condition(double y)		{ return cos(y);			} // функция, задающая граничное условие при x = 0
	double right_condition(double y)	{ return cos(y);			} // функция, задающая граничное условие при x = X
	double bottom_condition(double x)	{ return sin(x) + 1;		} // функция, задающая граничное условие при y = 0
	double top_condition(double x)		{ return sin(x) - 1;		} // функция, задающая граничное условие при y = Y
	double f(double x, double y)		{ return -sin(x) - cos(y);	} // функция, задающая внешнее воздействие
};
#endif	// HEAD_TASK_2

//#endif // SOFT_GRADER










//
//
//
//void heat_dirichlet_sor(heat_task task, double* v)
//{
//	const int n = task.n;
//	const int m = task.m;
//	const double h_x = task.X / n;
//	const double h_y = task.Y / m;
//
//	//const double W = 1.5 + std::sin(PI * h_x / 2) / 2;
//	//const double W = 1;
//	const double W = 2 / (1 + std::sin(M_PI * std::min(h_x, h_y) / 2));
//
//	const double REV_X = 1 / h_x / h_x;
//	const double REV_Y = 1 / h_y / h_y;
//	const double D = 2 * (REV_X + REV_Y);
//
//
//	std::vector<std::vector<double>> result(n + 1, std::vector<double>(m + 1, 0.0));
//#pragma omp parallel for if(PARALLEL)
//	for (int i = 0; i <= n; ++i)
//	{
//		result[i][0] = task.bottom_condition(i * h_x);
//		result[i][m] = task.top_condition(i * h_x);
//	}
//#pragma omp parallel for if(PARALLEL)
//	for (int j = 0; j <= m; ++j)
//	{
//		result[0][j] = task.left_condition(j * h_y);
//		result[n][j] = task.right_condition(j * h_y);
//	}
//
//
//	std::vector<std::vector<double>> func(n + 1, std::vector<double>(m + 1, 0.0));
//#pragma omp parallel for if(PARALLEL)
//	for (int i = 0; i <= n; ++i)
//	{
//		for (int j = 0; j <= m; ++j)
//		{
//			func[i][j] = -task.f(i * h_x, j * h_y);
//		}
//	}
//
//
//	/*std::vector<std::vector<double>> b(n + 1, std::vector<double>(m + 1, 0.0));
//
//
//	for (int i = 0; i <= n; ++i)
//	{
//	for (int j = 0; j <= m; ++j)
//	{
//	b[i][j] = task.f(i * h_x, j * h_y);
//	}
//	}
//	for (int i = 1; i < n; ++i)
//	{
//	b[i][1] += REV_X * task.bottom_condition(i * h_x);
//	b[i][m - 1] += REV_X * task.top_condition(i * h_x);
//	}
//	for (int j = 1; j < m; ++j)
//	{
//	b[0][j] += REV_Y * task.left_condition(j * h_y);
//	b[n - 1][j] += REV_Y * task.right_condition(j * h_y);
//	}*/
//
//
//	auto new_result(result);
//	double eps = 1e-4;
//	double diff = eps;
//	int iter = 0;
//
//	/*while (diff >= eps)
//	{
//	diff = 0;
//	++iter;
//	for (int i = 1; i < n; ++i)
//	{
//	for (int j = 1; j < m; ++j)
//	{
//	new_result[i][j] = func[i][j] + REV_X * (result[i - 1][j] + result[i + 1][j]) +
//	REV_Y * (result[i][j - 1] + result[i][j + 1]);
//	new_result[i][j] /= -D;
//	diff = std::max(std::abs(new_result[i][j] - result[i][j]), diff);
//	}
//	}
//	std::swap(result, new_result);
//	}*/
//
//	//int max_iter = 10 * std::max(1, static_cast<int>(1 / h_x));
//	int max_iter = 2 / std::min(h_x, h_y) / M_PI * log(1 / 1e-7);
//	//std::cout << max_iter << std::endl;
//
//	for (int iter = 0; iter < max_iter; ++iter)
//	{
//		for (int k = 2; k < n + m - 1; ++k)
//		{
//#pragma omp parallel for if(PARALLEL)
//			for (int i = std::max(1, k - m + 1); i < std::min(k, n); ++i)
//			{
//				int j = k - i;
//				/*if (iter == 1)
//				{
//				std::cout << i << ' ' << j << std::endl;
//				}*/
//
//				result[i][j] = W * (REV_X * (result[i - 1][j] + result[i + 1][j]) +
//					REV_Y * (result[i][j - 1] + result[i][j + 1]) + func[i][j]) +
//					(1 - W) * D * result[i][j];
//				result[i][j] /= D;
//			}
//		}
//		/*for (int i = 1; i < n; ++i)
//		{
//		for (int j = 1; j < m; ++j)
//		{
//		new_result[i][j] = W * (func[i][j] + REV_X * (result[i - 1][j] + result[i + 1][j]) +
//		REV_Y * (result[i][j - 1] + result[i][j + 1]));
//		//new_result[i][j] += (1 - W) * D * result[i][j];
//		new_result[i][j] /= D;
//		}
//		}*/
//		//std::swap(result, new_result);
//	}
//
//	std::cout << max_iter << std::endl;
//
//#pragma omp parallel for if(PARALLEL)
//	for (int i = 0; i <= n; ++i)
//		for (int j = 0; j <= m; ++j)
//			v[i * (m + 1) + j] = result[i][j];
//
//	/*for (int j = m; j >= 0; --j)
//	{
//	for (int i = 0; i <= n; ++i)
//	{
//	std::cout << std::setw(13) << result[i][j] << " ";
//	//std::cout << std::setw(15) << b[i][j] << " ";
//	}
//	std::cout << std::endl;
//	}*/
//	std::cout << result[n / 2][m / 3] << std::endl;
//}


#include <omp.h>
#include <vector>
#include <algorithm>

#define PARALLEL 1
#define M_PI	3.14159265358979323846  

typedef std::vector<std::vector<double>> grid;

void heat_dirichlet_sor(heat_task task, double* v) {
	double h_x = task.X / task.n;
	double h_y = task.Y / task.m;

	double W = 2 / (1 + std::sin(M_PI * std::min(h_x, h_y) / 2));

	double REV_X = 1 / ( h_x * h_x );
	double REV_Y = 1 / ( h_y * h_y );
	double D = (REV_X + REV_Y) * 2;

	grid G(task.n + 1, std::vector<double>(task.m + 1, 0.0));
#pragma omp parallel for shared(task, h_x, G) if(PARALLEL)
	for (int i = 0; i <= task.n; ++i) {
		G[i][0] = task.bottom_condition(i * h_x);
		G[i][task.m] = task.top_condition(i * h_x);
	}
#pragma omp parallel for shared(task, h_y, G) if(PARALLEL)
	for (int j = 0; j <= task.m; ++j) {
		G[0][j] = task.left_condition(j * h_y);
		G[task.n][j] = task.right_condition(j * h_y);
	}

	grid func(task.n + 1, std::vector<double>(task.m + 1, 0.0));
#pragma omp parallel for shared(h_x, h_y, func, task) if(PARALLEL)
	for (int i = 0; i <= task.n; ++i) {
		for (int j = 0; j <= task.m; ++j)
			func[i][j] = -task.f(i * h_x, j * h_y);
	}

	int max_iter = (int)(2 / std::min(h_x, h_y) / M_PI * log(1 / 1e-7));

	for (int iter = 0; iter < max_iter; ++iter) {
		for (int k = 2; k < task.n + task.m - 1; ++k) {
#pragma omp parallel for shared(G, W, D, REV_X, REV_Y, func) if(PARALLEL)
			for (int i = std::max(1, k - task.m + 1); i < std::min(k, task.n); ++i) {
				int j = k - i;

				G[i][j] = W * (REV_X * (G[i - 1][j] + G[i + 1][j]) +
					REV_Y * (G[i][j - 1] + G[i][j + 1]) + func[i][j]) +
					(1 - W) * D * G[i][j];
				G[i][j] /= D;
			}
		}
	}

#pragma omp parallel for shared(v, G) if(PARALLEL)
	for (int i = 0; i <= task.n; ++i)
		for (int j = 0; j <= task.m; ++j)
			v[i * (task.m + 1) + j] = G[i][j];
}

int main(int argc, char **argv) {
#if HEAD_TASK_1
	heat_task task(2.0, 3.0, 256, 256);
#endif // HEAD_TASK_1

#if HEAD_TASK_2
	heat_task task(M_PI, M_PI, 64, 64);
#endif // HEAD_TASK_2

	double* v = new double[(task.m + 1) * (task.n + 1)];
	heat_dirichlet_sor(task, v);
	delete []v;

	return 0;
}