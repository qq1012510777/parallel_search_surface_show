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

class surface_show_size
{
public:
    surface_show_size(vector<Matrixfl> Matrix3R_distance, vector<Matrixfl> &Matrix3R_show_surface, size_t Nproc, double max_R_);
};

surface_show_size::surface_show_size(vector<Matrixfl> Matrix3R_distance, vector<Matrixfl> &Matrix3R_show_surface, size_t Nproc, double max_R_)
{
    cout << "-----start surface_show_size-----\n";
    size_t NZ = Matrix3R_distance.size();
    size_t NY = Matrix3R_distance[0].rows();
    size_t NX = Matrix3R_distance[0].cols();

    size_t arraysize = NX * NY * NZ;
    //bool tag10 = false, tag20 = false, tag40 = false, tag60 = false, tag80 = false, tag95 = false;

    for (size_t global = 0; global < arraysize; ++global)
    {
        size_t z = global % NZ;
        size_t y = (global / NZ) % NY;
        size_t x = global / (NY * NZ);

        if (Matrix3R_distance[z](y, x) != 0 && Matrix3R_distance[z](y, x) != 1)
        {
            //----------now alter the values which are in the sphere with R = Matrix3R_distance[z](y, x), centering at (z, y, z)-------
            //------ if that value at (z_1, y_1, x_1) is 1 (air) or 0 (bone), ignore it

            //-------if that value = -1 or < Matrix3R_distance[z](y, x), alter Matrix3R_distance[z_1](y_1, x_1) to Matrix3R_distance[z](y, x)

            double R = Matrix3R_distance[z](y, x) - 1;

            if (R + 1 > max_R_)
                R = max_R_ - 1;

            if (Matrix3R_show_surface[z](y, x) == -1) // unchecked pore pixel
            {
                Matrix3R_show_surface[z](y, x) = R + 1; //Matrix3R_distance[z](y, x);
            }
            else if (Matrix3R_show_surface[z](y, x) >= Matrix3R_distance[z](y, x))
            {
                // do nothing
            }
            else if (Matrix3R_show_surface[z](y, x) < Matrix3R_distance[z](y, x))
            {
                Matrix3R_show_surface[z](y, x) = R + 1;
            }

            int min_x_ = x - R;
            if (min_x_ < 0)
                min_x_ = 0;

            int max_x_ = x + R;
            if (max_x_ > (int)NX - 1)
                max_x_ = NX - 1;

            int min_y_ = y - R;
            if (min_y_ < 0)
                min_y_ = 0;

            int max_y_ = y + R;
            if (max_y_ > (int)NY - 1)
                max_y_ = NY - 1;

            int min_z_ = z - R;
            if (min_z_ < 0)
                min_z_ = 0;

            int max_z_ = z + R;
            if (max_z_ > (int)NZ - 1)
                max_z_ = NZ - 1;

                //cout << min_x_ << ", " << max_x_ << endl;
                //cout << min_y_ << ", " << max_y_ << endl;
                //cout << min_z_ << ", " << max_z_ << endl;
                // cubic
                // sphere in the cubic

                // first, find the maximum value

#pragma omp parallel for schedule(dynamic) num_threads(Nproc)
            for (int z_1 = min_z_; z_1 <= max_z_; ++z_1)
            {
                //cout << "the first loop: " << z_1 << endl;
                for (int y_1 = min_y_; y_1 <= max_y_; ++y_1)
                {

                    for (int x_1 = min_x_; x_1 <= max_x_; ++x_1)
                    {
                        //cout << x_1 << ", " << y_1 << ", " << z_1 << endl;
                        if (x_1 != (int)x || y_1 != (int)y || z_1 != (int)z) // do not check itself
                        {
                            //cout << x_1 << ", " << y_1 << ", " << z_1 << endl;
                            Vector3d A, B, C;
                            A << x, y, z;
                            B << x_1, y_1, z_1;
                            C = A - B;

                            if (C.norm() <= R) // if the point is within sphere
                            {
                                if (Matrix3R_show_surface[z_1](y_1, x_1) < R + 1 && Matrix3R_distance[z_1](y_1, x_1) > 1)
                                {
                                    Matrix3R_show_surface[z_1](y_1, x_1) = R + 1;

                                    // if (Matrix3R_distance[z_1](y_1, x_1) <= 1)
                                    // {
                                    //     cout << "change an air/bone: " << x_1 << ", ";
                                    //     cout << y_1 << ", " << z_1 << endl;
                                    // }
                                }
                            }
                        }
                    }
                }
            }
            //---------------
        }
        else if (Matrix3R_distance[z](y, x) == 0)
        {
            Matrix3R_show_surface[z](y, x) = 0;
        }
        else if (Matrix3R_distance[z](y, x) == 1)
        {
            Matrix3R_show_surface[z](y, x) = 1;
        }
        else
        {
            cout << "undefined behavior!\n";
            exit(0);
        }

        double finishing = (global / ((1.0) * arraysize)) * 100;

        cout << "-----finish surface_show_size----";
        printf("%0.2f", finishing);
        cout << "% -----\n";
    }
    cout << "-----finish surface_show_size-----\n";
};
