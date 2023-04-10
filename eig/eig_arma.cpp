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

using namespace arma;

tuple<vector<vector<double>>, vector<vector<double>>, vector<double>, vector<double>> eigen(vector<vector<double>> arr, const int NG) {

    vector<vector <double>> V(NG , vector<double> (NG)); 
    vector<vector <double>> VI(NG , vector<double> (NG)); 
    vector<double> D(NG);
    vector<double> DI(NG);

    vector<vector <double>> VR(NG , vector<double> (NG)); 
    vector<vector <double>> VII(NG , vector<double> (NG)); 
    vector<double> WR(NG);
    vector<double> WI(NG);

    cx_mat A(NG, NG);

    // A[0] = 3;
    // A[1] = -2;
    // A[2] = 4;
    // A[3] = -1;

    // sp_mat A(NG, NG);

    cx_vec DWR;
    cx_mat DVR;

    for (int i = 0; i < NG; i++){
        for (int j = 0; j < NG; j++) {
            A[i + NG * j] = arr[i][j];
        }
    }

    int index = 0;
    printf("[");
    for(int i = 0 ; i < NG ; i++){
        printf("[");
        for(int j = 0 ; j < NG-1 ; j++){
            printf("%.16e, ", A[index]);
            index++;
        }
        printf("%.6e],\n", A[index]);
        index++;
    }
    printf("]\n");

    eig_gen(DWR, DVR, A, "balance");
    // eigs_gen(eigval, eigvec, A, 0);  // find 5 eigenvalues/eigenvectors

    // DWR.print("WR");
    // eigvec.print("eigvec");

    // printf("wr[0] = %.6e\n", WR[0]);

    for(int i = 0 ; i < NG ; i++){
        WR[i] = DWR[i].real();

        for(int j = 0 ; j < NG ; j++){
            VR[j][i] = DVR[j + i*NG].real();
        }

        if(fabs(DWR[i].imag()) > pow(10, -20)){
            WI[i] = DWR[i].imag();

            for(int j = 0 ; j < NG ; j++){
                VII[j][i] = DVR[j + i*NG].imag();
            }
        }
    }

    // for(int i = 0 ; i < NG ; i++){
    //     D[i] = DWR[i].real();
    //     DI[i] = DWR[i].imag();

    //     for(int j = 0 ; j < NG ; j++){
    //         V[i][j] = DVR[i + j*NG].real();
    //         VI[i][j] = DVR[i + i*NG].imag();
    //     }
    // }

    // printArray(D);

    // for(int i = 0 ; i < NG ; i++){
    //     for(int j = 0 ; j < NG ; j++){
    //         printf("%.6e\t", VR[i][j]);
    //     }
    //     printf("\n");
    // }

    // printf("\n");

    // for(int i = 0 ; i < NG ; i++){
    //     for(int j = 0 ; j < NG ; j++){
    //         printf("%.6e\t", VL[i][j]);
    //     }   
    //     printf("\n");
    // }
    // printf("\n");

    int c = 0, c_index = 0, r_index = 0; 

    // for(int i = 0 ; i < NG ; i++){    
    //     if(WI[i] != 0) c++;
    // }

    // // for(int i = 0 ; i < NG ; i++){
    // //     for(int j = 0 ; j < NG ; j++){
    // //         printf("%.6e ", VR[i][j]);
    // //     }
    // //     printf("\n");
    // // }
    // // printf("\n");

    // // for(int i = 0 ; i < NG ; i++){
    // //     printf("\t\t D %d ---- \t %.6e\t +\t %.6e\t i\n", i+1, WR[i], WI[i]);
    // // }

    // for(int i = 0 ; i < NG ; i++){
    //     // if(WR[i].real() > 0 ) {
    //         D[r_index] = WR[i].real();

    //         for(int j = 0 ; j < NG ; j++){
    //             V[j][r_index] = VR[j + i*NG].real();
    //         }

    //         r_index++;
        
    // }

    // for(int i = 0 ; i < NG ; i++){
    //     if(WR[i].real() < 0) {
    //         D[r_index] = WR[i].real();

    //         for(int j = 0 ; j < NG ; j++){
    //             V[j][r_index] = VR[i + j*NG].real();
    //         }

    //         r_index++;
    //     }
    // }

    ///////////////////////////////////////////////////////////////////////////////////////////////////

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] > 0 && WI[i] == 0) {
            D[r_index] = WR[i];
            for(int j = 0 ; j < NG ; j++){
                V[j][r_index] = VR[j][i];            
            }
            r_index++;
        }
    }

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] > 0 && WI[i] > 0) {
            D[r_index] = WR[i];
            DI[r_index] = WI[i];

            for(int j = 0 ; j < NG ; j++){
                V[j][r_index] = VR[j][i];
                VI[j][r_index] = VII[j][i];


                V[j][r_index+1] = VR[j][i];
                VI[j][r_index+1] = -VR[j][i];
            }

            r_index++;
        }
    }

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] > 0 && WI[i] < 0) {
            D[r_index] = WR[i];
            DI[r_index] = WI[i];

            for(int j = 0 ; j < NG ; j++){
                V[j][r_index] = VR[j][i];
                VI[j][r_index] = VII[j][i];


                V[j][r_index + 1] = VR[j][i];
                VI[j][r_index + 1] = -VII[j][i];
            }

            r_index++;
        }
    }

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] < 0 && WI[i] == 0) {
            D[r_index] = WR[i];

            for(int j = 0 ; j < NG ; j++){
                V[j][r_index] = VR[j][i];
            }

            r_index++;
        }
    }

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] < 0 && WI[i] > 0) {
            D[r_index] = WR[i];
            DI[r_index] = WI[i];

            for(int j = 0 ; j < NG ; j++){
                V[j][r_index] = VR[j][i];
                VI[j][r_index] = VII[j][i];


                V[j][r_index+1] = VR[j][i];
                VI[j][r_index+1] = -VR[j][i];
            }

            r_index++;
        }
    }

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] < 0 && WI[i] < 0) {
            D[r_index] = WR[i];
            DI[r_index] = WI[i];

            for(int j = 0 ; j < NG ; j++){
                V[j][r_index] = VR[j][i];
                VI[j][r_index] = VII[j][i];


                V[j][r_index+1] = VR[j][i];
                VI[j][r_index+1] = -VR[j][i];
            }

            r_index++;
        }
    }

    // free(VV);
    
    return {V, VI, D, DI};
}
