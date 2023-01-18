#pragma once

#include "load_3rank_data.h"
#include <Eigen/Dense>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

class cubic_searching
{
public:
    vector<MatrixXs> Matrix_distance; // [z](y , x)
public:
    cubic_searching(const vector<MatrixXs> Matrix3R, size_t maxE_1, size_t bone, size_t pore, size_t air, vector<Matrixfl> &Matrix3R_distance, size_t Nproc);
};

cubic_searching::cubic_searching(const vector<MatrixXs> Matrix3R, size_t maxE_1, size_t bone, size_t pore, size_t air, vector<Matrixfl> &Matrix3R_distance, size_t Nproc)
{
    //-------------------------------------------------------------------------------------
    //a new parallization scheme
    //-------------------------------------------------------------------------------------
    cout << "-----start serarching-----\n";
    size_t NZ = Matrix3R.size();
    size_t NY = Matrix3R[0].rows();
    size_t NX = Matrix3R[0].cols();
    double yu = 0;
    size_t arraysize = NX * NY * NZ;
    //bool tag10 = false, tag20 = false, tag40 = false, tag60 = false, tag80 = false, tag95 = false;

#pragma omp parallel for schedule(dynamic) num_threads(Nproc)
    for (size_t global = 0; global < arraysize; ++global)
    {
        size_t z = global % NZ;
        size_t y = (global / NZ) % NY;
        size_t x = global / (NY * NZ);

        //size_t x = 388;
        //size_t y = 238;
        //size_t z = 1;

        //cout << Matrix3R[z](y, x)  << endl;

        if (Matrix3R[z](y, x) == bone) // 0
        {
            Matrix3R_distance[z](y, x) = bone;
        }
        else if (Matrix3R[z](y, x) == air) // 1
        {
            Matrix3R_distance[z](y, x) = pore; //1
        }
        else if (Matrix3R[z](y, x) == pore) // 2
        {
            size_t last_minX = x;
            size_t last_maxX = x;
            size_t last_minY = y;
            size_t last_maxY = y;
            size_t last_minZ = z;
            size_t last_maxZ = z;

            double dist_pore_2_bone = 1e10;
            bool ui = false;

            size_t maxE = maxE_1;
            // cubic searching starts-----------------------------------
            for (size_t i = 1; i <= maxE; ++i)
            {
                int minX = x - i;
                int maxX = x + i;

                int minY = y - i;
                int maxY = y + i;

                int minZ = z - i;
                int maxZ = z + i;

                if (minX < 0)
                    minX = 0;

                if (maxX > (int)NX - 1)
                    maxX = NX - 1;

                if (minY < 0)
                    minY = 0;

                if (maxY > (int)NY - 1)
                    maxY = NY - 1;

                if (minZ < 0)
                    minZ = 0;

                if (maxZ > (int)NZ - 1)
                    maxZ = NZ - 1;

                for (int z_1 = minZ; z_1 <= maxZ; ++z_1)
                {
                    for (int y_1 = minY; y_1 <= maxY; ++y_1)
                    {
                        for (int x_1 = minX; x_1 <= maxX; ++x_1)
                        {

                            if ((z_1 > (int)last_maxZ || z_1 < (int)last_minZ) ||
                                (y_1 > (int)last_maxY || y_1 < (int)last_minY) ||
                                (x_1 > (int)last_maxX || x_1 < (int)last_minX))
                            {
                                /*
                                if (z == 0)
                                    cout << x_1 << ", " << y_1 << ", " << z_1 << ": " << Matrix3R[z_1](y_1, x_1) << endl;
                                */
                                if (Matrix3R[z_1](y_1, x_1) == bone) // 1
                                {
                                    double f = (Vector3d{(double)x, (double)y, (double)z} - Vector3d{(double)x_1, (double)y_1, (double)z_1}).norm();

                                    //cout << "1g f: " << f << ", dist_pore_2_bone: " << dist_pore_2_bone << endl;
                                    if (f < dist_pore_2_bone)
                                    {
                                        //cout << "2g f: " << f << ", dist_pore_2_bone: " << dist_pore_2_bone << endl;
                                        dist_pore_2_bone = f;

                                        Matrix3R_distance[z](y, x) = dist_pore_2_bone + 1;
                                        //cout << x_1 << ", " << y_1 << ", " << z_1 << ": " << Matrix3R[z_1](y_1, x_1) << endl;
                                        //cout << "\tabs(dist_pore_2_bone - i) " << abs(dist_pore_2_bone - i) << endl;
                                        if (ui == false && abs(dist_pore_2_bone - i) < 1e-1) // distance is the radius of the sphere
                                        {
                                            //cout << 2 << endl;
                                            goto here1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (2 == 1)
                {
                here1:;
                    Matrix3R_distance[z](y, x) = dist_pore_2_bone + 1;
                    break;
                }

                if (ui == true)
                    break; // break cubic searching of this point

                last_minX = minX;
                last_maxX = maxX;
                last_minY = minY;
                last_maxY = maxY;
                last_minZ = minZ;
                last_maxZ = maxZ;

                if (i == maxE && dist_pore_2_bone == 1e10)
                {
                    Matrix3R_distance[z](y, x) = maxE + 1;
                }

                if ((dist_pore_2_bone != 1e10) && ((dist_pore_2_bone - i > 0.1))) // found bond may be not within the sphere (R = i)
                {

                    ui = true;
                    size_t maxR_o = ceil(dist_pore_2_bone - i) + i;

                    if (maxR_o <= maxE)
                    {
                        maxE = maxR_o;
                        i = maxE - 1;
                    }
                    else
                    {
                        i = maxR_o - 1;
                        maxE = maxR_o;
                    }
                    //cout << i << ", " << maxE << endl;
                }

                if (dist_pore_2_bone - i < -0.01)
                {
                    cout << "error\n";
                    exit(0);
                }
            }

            
            //cout << Matrix3R_distance[z](y, x) << endl;

            // cubic searching finished-----------------------------------
        }
        else
        {
            cout << "undefine material!!!\n";
            exit(0);
        }

        //if (1)

        /*
#pragma omp critical
        {

            yu++;

            double finishing = (yu / ((1.0) * arraysize)) * 100;
            if (1)
            {
                cout << "-----finish ";
                printf("%0.2f", finishing);
                cout << "% -----\n";
            }
        }*/
    }
    cout << "-----finish serarching-----\n";
};
