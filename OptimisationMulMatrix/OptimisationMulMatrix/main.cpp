#include "windows.h"
#include "MatrixAlgorithms.h"

#define SIZE_MATRIX_SHOW    5
#define RANDOM_RANGE        5.0

#define FILE_RESULT_LAUNCH  "ÎptResult"

int main(int argc, char** argv) {
  double start = 0;
  int count_it = 0;

  std::vector<int> arr_matrix_size = { 100, 200, 500, 1000, 1500, 2000 };
  std::vector<double> run_time_before_opt, run_time_after_opt;

  MSG("\n========== Matrix Multiplication Algorithm =========");
  for (auto it : arr_matrix_size) {
    printf("=============== MATRIX SIZE: %dx%d ===============\n\n", it, it);

    double *left_matrix = init_matrix(it);
    double *rigth_matrix = init_matrix(it);
    double *res_matrix = init_matrix(it);

    fill_matrix(left_matrix, it, RANDOM_RANGE);
    fill_matrix(rigth_matrix, it, RANDOM_RANGE);

    start = omp_get_wtime();
    {
      multiplication_matrix(res_matrix, left_matrix, rigth_matrix, it);
    }
    run_time_before_opt.push_back(omp_get_wtime() - start); MSG("Not optimization");
    check_result_mul_matrix(res_matrix, left_matrix, rigth_matrix, it, SIZE_MATRIX_SHOW);
    printf("Run Time: %.5f\n\n", run_time_before_opt[count_it]);

    start = omp_get_wtime();
    {
      optimal_multiplication_matrix(res_matrix, left_matrix, rigth_matrix, it);
    }
    run_time_after_opt.push_back(omp_get_wtime() - start); MSG("Optimization");
    check_result_mul_matrix(res_matrix, left_matrix, rigth_matrix, it, SIZE_MATRIX_SHOW);
    printf("Run Time: %.5f\n\n", run_time_after_opt[count_it]);

    printf("Speedup: %.5f\n", run_time_before_opt[count_it] / run_time_after_opt[count_it]);
    ++count_it;

    DESTROY_MATRIX(left_matrix);
    DESTROY_MATRIX(rigth_matrix);
    DESTROY_MATRIX(res_matrix);
  }

  print_result_experiment_to_file(arr_matrix_size, run_time_before_opt, run_time_after_opt, "ÎptResult");

  MSG("\n\nTo end, press any key");
  _getch();

  return EXIT_SUCCESS;
}
