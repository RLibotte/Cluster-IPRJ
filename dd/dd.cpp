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
#include <tuple>
#include <omp.h>
#include <time.h>
#include <vector>
#include <string>
#include <fstream>

#include "./create_arr.h"

#include "../print.cpp"
#include "../get_GQ.cpp"
#include "../get_dev.cpp"
#include "../get_phi.cpp"
#include "../get_reflexive_bounds.cpp"
#include "../print_output.cpp"

#include "./create_array.cpp"
#include "./get_psi_m.cpp"
#include "./get_SS.cpp"
#include "./sweep.cpp"
#include "./nodes_array.cpp"

#include "../models/model_1.cpp"

void dd(const int N, double error){

    //////////////////////////////////////////////////////////////////////////////////////////

    get_GQ(N, mi, w);
    get_model(R, Z, G, L, h_R, q_R, sigmat_R, sigmas0_R, sigmas1_R, sigmat_z, sigmas0_z, sigmas1_z, nodes, zone_config, nxf, lb, rb);
    create_arrays(N, G, nxf, L, lb, rb, psi, psi_m, phi, phi_old, SS);
    nodes_arr(N, L, R, G, nxf, nodes, sigmat_R, sigmas0_R, sigmas1_R, q_R, h_R, sigmat, q, sigmas0, sigmas1, h, psi_m, SS, pos_arr);

    /////////////////////////////////////////////////////////////////////////////////////////////

    do{

        psi_m = get_psi_m(psi, psi_m, N*G, nxf);
        SS = get_SS(SS, sigmas0, sigmas1, mi, w, psi_m, nxf, G, N, L);

        psi = sweep_lr(psi, mi, sigmat, SS, q, h, nxf, G, N);
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        psi_m = get_psi_m(psi, psi_m, N*G, nxf);
        SS = get_SS(SS, sigmas0, sigmas1, mi, w, psi_m, nxf, G, N, L);
        
        psi = sweep_rl(psi, mi, sigmat, SS, q, h, nxf, G, N);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        if(lb == (-1)) psi = reflexive_l(psi, N, G);
        if(rb == (-1)) psi = reflexive_r(psi, N, G, nxf);

        phi = get_phi(phi, psi, w, nxf, G, N);
        dev = get_dev(phi, phi_old, G, nxf);

        phi_old = phi;

        a++;

        printf("%d --- %.6e \n", a, dev);

    // } while (a < 10);
    } while (dev > error);

    int index = 0;
    printf("\n--- Phi ---\n");
    for(int j = 0 ; j < G ; j++){
        for(int i = 0 ; i < nxf+1 ; i++){
            printf("%.6f ", phi[j][i]);
        }
        printf("\n");
    }
    for(int j = 0 ; j < G ; j++){
        for(int i = 0 ; i < R ; i++){
            printf("%.6e ", phi[j][index]);
            index += nodes[i]; 
        }
        printf("%.6e \n\n", phi[j][index]);
        index = 0;
    }

        // for(int i = 0 ; i < R ; i++){
        //     printf("%.6e & ", phi[0][index]);
        //     index += nodes[i]; 
        // }
        // printf("%.6e  \\\\ \n\n", phi[0][index]);
        // index = 0;

        // for(int i = 0 ; i < R ; i++){
        //     printf("%.6e & ", phi[G-1][index]);
        //     index += nodes[i]; 
        // }
        // printf("%.6e \\\\ \n\n", phi[G-1][index]);
        // index = 0;


    pos_arr = get_pos_arr(R, nxf, h, nodes);
    print_output(G, nxf, phi, pos_arr);
    
}


int main (){

    clock_t tStart = clock();

    const int num_threads = 1;
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        printf("%d THREADS\n", omp_get_num_threads());
    }

    const int N = 2096;
    double error = pow(10, -6);

    dd(N, error);

    printf("Time taken: %.6f s\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}