////////////////////////////////////////////////////////////////
//                                                            //
//         Universidade do Estado do Rio de Janeiro           //
//                    Instituto Politécnico                   //
//                           PPGMC                            //
//                                                            //
//              Simulador Transporte de Nêutrons              //
//                                                            //
//               Autor: Rafael Barbosa Libotte                //
//                                                            //
////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <vector>
#include <fstream>

#include "./create_arr.h"

#include "../eig/Alglib/linalg.h"
#include "../eig/Alglib/solvers.h"
#include "../eig/Alglib/ap.h"

#include "../coarse_mesh/create_array.cpp"
#include "../coarse_mesh/get_arrays.cpp"
// #include "../coarse_mesh/get_psi_m.cpp"
// #include "../coarse_mesh/get_SS.cpp"


// #include "../print.cpp"
#include "../eig/eig_alg.cpp"
#include "../get_eig.cpp"
#include "../get_GQ.cpp"
#include "../get_dev.cpp"
#include "../get_phi.cpp"
#include "../get_inverse.cpp"
#include "../get_alpha.cpp"
#include "../get_reflexive_bounds.cpp"
// #include "../get_std_dev.cpp"
#include "../print_output.cpp"
// #include "../print_ref.cpp"

#include "./get_psi.cpp"

#include "../models/model_1.cpp"

void sdm(const int N, double error){

    get_GQ(N, mi, w);
    get_model(R, Z, G, L, h_R, q_R, sigmat_R, sigmas0_R, sigmas1_R, sigmat_z, sigmas0_z, sigmas1_z, nodes, zone_config, nxf, lb, rb);
    create_arrays(N, G, R, L, lb, rb, psi, psi_m, phi, phi_old, SS, alpha, alpha_temp, avg, D_eig, DI_eig, V_eig, VI_eig, ni, ni_i, V, VI);
    get_eig(R, Z, G, N, sigmat_z, sigmas0_z, sigmas1_z, mi, w, zone_config, ni, ni_i, V, VI);

    const int NG = N*G;

    for(int i = 0 ; i < R ; i++){
        for(int j = 0 ; j < NG ; j++){
            ni[j][i] = 1/ni[j][i];
        }
    }

    get_arrays(N, R, G, V, ni, h_R, A_p, q_R, q_p, psi_p_temp, temp, A_p_temp, psi_p, Am);
    get_inverse(NG, R, Am, Am_inv);

    // ///////////////////////////////////////////////////////////////////////////////////////////

    vector<vector<double>> Am_temp (NG, vector<double> (N*G));
    vector<double> psi_temp(NG);

    for(int i = 0 ; i < iter ; i++){
        a = 0;
        dev = 0;
        // clock_t tStart = clock();

        // for(int l = 0 ; l < R ; l++){
        //     for(int j = 0 ; j < G ; j++){
        //         for(int k = 0 ; k < N/2 ; k++){
        //             psi[k + j*N][l+1] = 0;
        //             psi[k+N/2 + j*N][l] = 0;
        //         }
        //     }
        // }

        // for(int j = 0 ; j < N/2 ; j++){
        //     psi[j][0] = 1;
        // }

        do{
            for(int k = 0 ; k < R ; k++){
                alpha[k] = get_alpha(N, G, k, psi, psi_p, Am_inv[k]);

                get_psi_l(N, G, NG, k, psi, ni, Am_temp, V, h_R, psi_p, alpha[k], psi_temp);
                get_psi_r(N, G, NG, k, psi, ni, Am_temp, V, h_R, psi_p, alpha[k], psi_temp);
            }

            if(lb == (-1)) psi = reflexive_l(psi, N, G);
            if(rb == (-1)) psi = reflexive_r(psi, N, G, R);

            phi = get_phi(phi, psi, w, R, G, N);
            dev = get_dev(phi, phi_old, G, R);
            phi_old = phi;

            a++;
            
            // printf("\npsi\n\n");
            // print2dArray(psi);

            // printf("\nphi \n\n");
            // print2dArray(phi);

            printf("%d \n", a);

        } while (a < 2);
        // } while (dev > error);

        // sum += (double)(clock() - tStart)/CLOCKS_PER_SEC;
        // avg[i] = (double)(clock() - tStart)/CLOCKS_PER_SEC;
        // printf("Iterative Process (%d): %.6f\n", i+1, (double)(clock() - tStart)/CLOCKS_PER_SEC);
    }

    // printf("\n\nAverage time - iterative process = %.6f s\n", sum/iter);
    // get_std_dev(iter, sum, avg);

    printf("\npsi\n\n");
    print2dArray(psi);

    printf("\nphi\n\n");
    print2dArray(phi);

    // printf("%d iterations\n", a);

    pos_arr = get_pos_arr_coarse(R, h_R);
    print_output(G, R, phi, pos_arr);

}


int main (){

    clock_t tStart = clock();

    const int num_threads = 1;
    omp_set_num_threads(num_threads);
    #pragma omp parallel 
    {
        printf("%d THREADS\n", omp_get_num_threads());
    }

    const int N = 4;
    double error = pow(10, -6);

    sdm(N, error);

    printf("Time taken: %.6f s\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}
