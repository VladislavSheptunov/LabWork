#include <omp.h>
#include <assert.h>
#include <vector>
#include <algorithm>

#define PARALLEL 1
#define SOFT_GRADER 0

#define M_PI	3.14159265358979323846

typedef std::vector< std::vector< double > > grid;

#if !SOFT_GRADER
#include <conio.h>
#include <fstream>
#define MSG(msg) printf("%s\n", ##msg)
// TEST TASK
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

static inline void check_convergence(const heat_task &task, double* result) {
	// @TODO: Make check
	MSG("SUCCESSFUL");
}

static inline void print_result_experiment_to_file(	const std::vector<int> &size_matrix, 
													const std::vector<int> &thread, 
													const std::vector<double> &speedup) 
{
	system("del ./TableSpeedUp.csv");

	std::ofstream TableSpeedUp;
	TableSpeedUp.open("TableSpeedUp.csv");
	//TableSpeedUp.open("TableEfficiency.csv");
	TableSpeedUp << "Count thread\\Size grid;";

	for (auto it : size_matrix)
		TableSpeedUp << it << ";";

	TableSpeedUp << "\n";

	for (int it = 0; it < thread.size(); it++) {
		TableSpeedUp << thread[it] << ";";
		for (int i = 0; i < size_matrix.size(); i++)
			TableSpeedUp << speedup[it + i * thread.size()] << ";";
		TableSpeedUp << "\n";
	}

	TableSpeedUp.close();
}

#endif // SOFT_GRADER

void heat_dirichlet_sor(heat_task task, double* v) {
	const double h_x = task.X / task.n;
	const double h_y = task.Y / task.m;

	double w_opt = 2 / (1 + std::sin(M_PI * std::min(h_x, h_y) / 2));

	double REV_X = 1 / ( h_x * h_x );
	double REV_Y = 1 / ( h_y * h_y );
	double alpha = (REV_X + REV_Y) * 2;

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
				//  Уравнения разностной схемы 
				G[i][j] = w_opt * (REV_X * (G[i - 1][j] + G[i + 1][j]) + REV_Y * (G[i][j - 1] + G[i][j + 1]) + func[i][j]) + (1 - w_opt) * alpha * G[i][j];
				G[i][j] /= alpha;
			}
		}
	}

#pragma omp parallel for shared(v, G) if(PARALLEL)
	for (int i = 0; i <= task.n; ++i)
		for (int j = 0; j <= task.m; ++j)
			v[i * (task.m + 1) + j] = G[i][j];
}

#if !SOFT_GRADER
int main(int argc, char **argv) {
	// Starting from the debug
	std::vector<int> arr_grid_size = { 100, 250, 500, 750, 1000, 1250, 1500 };
	std::vector<int> count_thread;
	std::vector<double> speedUp;
	int max_threads = 4; // omp_get_max_threads();
	for (int i = 1; i < max_threads + 1; i++)
		count_thread.push_back(i);

	double grid_x = 2.0, // ширина пластины
		   grid_y = 3.0; // длина пластины
	
	double *start_time = new double[max_threads + 1];
	double *diff_time = new double[max_threads + 1];
	double *speedup = new double[max_threads + 1];
	double *efficiency = new double[max_threads + 1];
	//--

	for (auto it : arr_grid_size) {
		printf("\n========== GRID SIZE: %dx%d PLATE SEZE %.1fx%.1f =========\n\n", it, it, grid_x, grid_y);
		omp_set_num_threads(1);
		heat_task task(grid_x, grid_y, it, it);
		double* v = new double[(task.m + 1) * (task.n + 1)];
		
		start_time[0] = omp_get_wtime();
		heat_dirichlet_sor(task, v);
		diff_time[0] = omp_get_wtime() - start_time[0];
		speedup[0] = 1;
		efficiency[0] = speedup[0] / 1;

		printf("The result of the sequential algorithm:\n");
		printf("Run time:   %.5f \n", diff_time[0]);
		printf("Speedup:    %.5f \n", speedup[0]);
		printf("Efficiency: %.5f \n\n", efficiency[0]);

		check_convergence(task, v);

		for (int i = 1; i < max_threads + 1; i++) {
			if (i == 1)
				continue;
			printf("\n========== GRID SIZE: %dx%d PLATE SEZE %.1fx%.1f =========\n\n", it, it, grid_x, grid_y);
			omp_set_num_threads(i);
			std::fill_n(v, (task.m + 1) * (task.n + 1), 0.0);
			start_time[i] = omp_get_wtime();
			heat_dirichlet_sor(task, v);
			diff_time[i] = omp_get_wtime() - start_time[i];
			speedup[i] = diff_time[0] / diff_time[i];
			efficiency[i] = speedup[i] / i;
			printf("\nThe result of the parallel algorithm (number of threads %d):\n", i);
			printf("Run time:   %.5f \n", diff_time[i]);
			printf("Speedup:    %.5f \n", speedup[i]);
			printf("Efficiency: %.5f \n\n\n", efficiency[i]);
			speedUp.push_back(speedup[i]);
			//speedUp.push_back(efficiency[i]);
			check_convergence(task, v);
		}

		delete[]v;
	}

	print_result_experiment_to_file(arr_grid_size, count_thread, speedUp);

	delete[]start_time;
	delete[]diff_time;
	delete[]speedup;
	delete[]efficiency;

	MSG("To end, press any key");
	_getch();

	return EXIT_SUCCESS;
}
#endif // SOFT_GRADER