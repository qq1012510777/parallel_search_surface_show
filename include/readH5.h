#pragma once
#include "H5Cpp.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include "hdf5.h"
#include "Eigen/Dense"
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace H5;

class readH5
{
public:
    vector<int> Data;
    int Rank;
    vector<size_t> Dim;

public:
    readH5(string FileName, string DatasetName);
};

readH5::readH5(string FileName, string DatasetName)
{
    cout << "-----start to load data-----\n";
    H5File h5_read_file(FileName, H5F_ACC_RDONLY);

    // access the required dataset by path name
    DataSet dset = h5_read_file.openDataSet(DatasetName);

    // get the dataspace
    DataSpace dspace = dset.getSpace();

    // get the dataset type class
    //H5T_class_t type_class = dset.getTypeClass();

    // get the size of the dataset
    int rank = dspace.getSimpleExtentNdims(); //
    cout << "Rank: " << rank << endl;
    Rank = rank;

    hsize_t dims[rank] = {};
    rank = dspace.getSimpleExtentDims(dims);

    Dim.resize(rank);
    copy(dims, dims + rank, Dim.begin());

    cout << "Datasize: "; // this is the correct number of values
    for (size_t i = 0; i < (size_t)rank; ++i)
        cout << dims[i] << ", ";
    cout << endl;
    DataSpace memspace(rank, dims);

    size_t arraysize = 1;
    for (size_t i = 0; i < (size_t)rank; ++i)
        arraysize = arraysize * dims[i];

    int *p1DArray = new int[arraysize]();
    cout << "arraysize: " << arraysize << endl;
    dset.read(p1DArray, PredType::NATIVE_INT, memspace, dspace);

    this->Data.resize(arraysize);
    copy(p1DArray, p1DArray + arraysize, Data.begin());

    delete[] p1DArray;
    p1DArray = NULL;
    h5_read_file.close();
    cout << "-----finish to load data-----\n";
};
