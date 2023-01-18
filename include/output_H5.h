#pragma once
#include "H5Cpp.h"
#include "cubic_searching.h"
#include "hdf5.h"
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace H5;

class output_H5
{
public:
    output_H5(string filename, string dataset_name, vector<Matrixfl> Matrix3R_distance, readH5 r);
};

output_H5::output_H5(string filename, string dataset_name, vector<Matrixfl> Matrix3R_distance, readH5 r)
{
    cout << "-----start outputing data-----\n";
    H5File file(filename, H5F_ACC_TRUNC);

    hsize_t dimsf[r.Rank];

    size_t arraysize = 1;
    for (size_t i = 0; i < (size_t)r.Rank; ++i)
    {
        dimsf[i] = r.Dim[i];
        arraysize = arraysize * r.Dim[i];
    }

    DataSpace dataspace(r.Rank, dimsf);

    IntType datatype(PredType::NATIVE_FLOAT);

    datatype.setOrder(H5T_ORDER_LE);

    DataSet dataset = file.createDataSet(dataset_name, datatype, dataspace);

    _Float32 *buffer = new _Float32[arraysize]();

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

                        buffer[globalIndex] = Matrix3R_distance[z](y, x);
                        /*
                        if (buffer[globalIndex] != 0)
                        {
                            if (buffer[globalIndex] - (int)buffer[globalIndex] > 0)
                            {
                                cout << buffer[globalIndex] - (int)buffer[globalIndex] << endl;
                            }
                        }
                        */
                        
                    }
                }
            }
        }
    }

    dataset.write(buffer, PredType::NATIVE_FLOAT);

    file.close();
    delete[] buffer;
    buffer = NULL;
    cout << "-----finish outputing data-----\n";
};