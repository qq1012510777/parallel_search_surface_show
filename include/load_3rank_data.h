#pragma once
#include "readH5.h"
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace H5;

typedef Matrix<size_t, Dynamic, Dynamic> MatrixXs;
typedef Matrix<float, Dynamic, Dynamic> Matrixfl;

class load_3rank_data
{
public:
    load_3rank_data(vector<MatrixXs> &Matrix3R, readH5 r);
};

load_3rank_data::load_3rank_data(vector<MatrixXs> &Matrix3R, readH5 r)
{
    size_t NT = 1, NZ = r.Dim[1], NY = r.Dim[2], NX = r.Dim[3], NC = 1;
    for (size_t t = 0; t < NT; ++t)
    {
        for (size_t z = 0; z < NZ; ++z)
        {
            for (size_t y = 0; y < NY; ++y)
            {
                for (size_t x = 0; x < NX; ++x)
                {
                    for (size_t c = 0; c < NC; ++c)
                    {
                        // t z y x c
                        size_t globalIndex = c + x * NC + y * NX * NC + z * NY * NX * NC + t * NZ * NY * NX * NC;

                        Matrix3R[z](y, x) = r.Data[globalIndex];
                    }
                }
            }
        }
    }
};