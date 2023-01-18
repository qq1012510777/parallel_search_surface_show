#include "../include/cubic_searching.h"
#include "../include/load_3rank_data.h"
#include "../include/output_H5.h"
#include "../include/readH5.h"
#include "../include/surface_show_size.h"
#include "H5Cpp.h"
#include "hdf5.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace H5;

typedef Matrix<size_t, Dynamic, Dynamic> MatrixXs;

typedef Matrix<float, Dynamic, Dynamic> Matrixfl;

int main()
{
    readH5 r{"../inp/cut1.h5", "/data"};
    // now we know the t , z, y , x, c

    MatrixXs m = MatrixXs::Zero(r.Dim[2], r.Dim[3]);
    vector<MatrixXs> Matrix3R{r.Dim[1], m};
    // Matrix3R's element_index is [z](y, x)
    load_3rank_data l{Matrix3R, r};
    r.Data.clear();

    size_t maxE = 1e10;

    Matrixfl ks = Matrixfl::Zero(r.Dim[2], r.Dim[3]);

    vector<Matrixfl> Matrix3R_distance{r.Dim[1], ks};

    Matrixfl kl = Matrixfl::Zero(r.Dim[2], r.Dim[3]);
    for (size_t i = 0; i < r.Dim[2]; ++i)
        for (size_t j = 0; j < r.Dim[3]; ++j)
            kl(i, j) = -1;

    vector<Matrixfl> Matrix3R_show_surface{r.Dim[1], kl};

    size_t NUM_proc = 10;
    size_t backbone = 0;
    size_t air = 2;
    size_t pore = 1;
    cubic_searching c{Matrix3R, maxE, backbone, pore, air, Matrix3R_distance, NUM_proc}; // NUM_proc is the number of CPU
    output_H5 output1{"distance_of_pore_to_bone.h5", "/data", Matrix3R_distance, r};

    double max_R = 20;
    surface_show_size sss{Matrix3R_distance, Matrix3R_show_surface, NUM_proc, max_R};

    output_H5 output2{"size_of_pore.h5", "/data", Matrix3R_show_surface, r};
    
    return 0;
};