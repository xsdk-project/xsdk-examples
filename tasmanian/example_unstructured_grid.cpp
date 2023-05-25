
#include "Tasmanian.hpp"
#include <random>

using namespace std;

void sparse_grids_example_11(){
    cout << "\n---------------------------------------------------------------------------------------------------\n";
    cout << std::scientific; cout.precision(4);
    cout << "Construction of a sparse grid surrogate but using unstructured data\n\n";

    int const num_inputs = 2;

    // using random points to test the error
    int const num_test_points = 1000;
    std::vector<double> test_points(num_test_points * num_inputs);
    std::minstd_rand park_miller(42);
    std::uniform_real_distribution<double> domain(-1.0, 1.0);
    for(auto &t : test_points) t = domain(park_miller);

    // computes the error between the gird surrogate model and the actual model
    // using the test points, finds the largest absolute error
    auto get_error = [&](TasGrid::TasmanianSparseGrid const &grid,
                         std::function<void(double const x[], double y[], size_t)> model)->
        double{
            std::vector<double> grid_result;
            grid.evaluateBatch(test_points, grid_result);
            double err = 0.0;
            for(int i=0; i<num_test_points; i++){
                double model_result; // using only one output
                model(&test_points[i*num_inputs], &model_result, 0);
                err = std::max(err, std::abs(grid_result[i] - model_result));
            }
            return err;
        };

    // using a simple model
    auto model = [](double const x[], double y[], size_t)->
        void{ y[0] = std::exp(-x[0]*x[0] -x[1]-x[1]); };

    auto grid = TasGrid::makeGlobalGrid(num_inputs, 1, 4, TasGrid::type_level, TasGrid::rule_clenshawcurtis);

    // generate random data for the inputs, and compute the corresponding outputs
    int const num_data_points = 2000;
    std::vector<double> data_input(num_inputs * num_data_points);
    std::vector<double> data_output(num_data_points);

    for(auto &d : data_input) d = domain(park_miller);
    for(int i=0; i<num_data_points; i++) model(&data_input[i * num_inputs], &data_output[i], 0);

    // check if capability is available
    if (not grid.isAccelerationAvailable(TasGrid::accel_cpu_blas) and
        not grid.isAccelerationAvailable(TasGrid::accel_gpu_cuda) and
        not grid.isAccelerationAvailable(TasGrid::accel_gpu_magma)){
        cout << "Skipping the example, BLAS, CUDA/ROCm, or MAGMA acceleration required.\n";
        return;
    }

    // accel_cpu_blas: works on the CPU and can utilize all availabel RAM
    // accel_gpu_cuda: works on the GPU but it is resticted to the case
    //                 where the data can fit in GPU memory
    // accel_gpu_magma: works out-of-core, the data is stored in CPU RAM
    //                  while computations are still done on the GPU
    grid.enableAcceleration(TasGrid::accel_gpu_magma);

    // constructs the grid, depending on the amoutn of data data,
    // the side of the grid and the enabled acceleration
    // this can take significant amount of time
    TasGrid::loadUnstructuredDataL2(data_input, data_output, 1.E-4, grid);

    cout << "Using construction from unstructured (random) data\n";
    cout << "   approximatino error = " << get_error(grid, model) << "\n\n";
}

int main(int, char**){

    sparse_grids_example_11();

    return 0;
}
