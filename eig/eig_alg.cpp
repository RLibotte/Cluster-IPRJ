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

void eigen(std::vector<std::vector<double>> arr, const int NG, 
    std::vector<double> &D_eig, std::vector<double> &DI_eig,
    std::vector<std::vector<double>> &V_eig, std::vector<std::vector<double>> &VI_eig){

    int VNeeded = 3;
    bool Conv;
    double *VV; 
    alglib::real_2d_array AA, VR, VL;
    alglib::real_1d_array WR, WI;

    VV = new double[NG * NG];
    VR.setlength(NG, NG);
    VL.setlength(NG, NG);
    WR.setlength(NG);
    WI.setlength(NG);

    int k = 0;
    for (int i = 0; i < NG; i++){
        for (int j = 0; j < NG; j++) {
            VV[j + NG * i] = arr[i][j];
        }
    }

    AA.setcontent(NG, NG, VV);
    Conv = alglib::rmatrixevd(AA, NG, VNeeded, WR, WI, VL, VR);

    int c = 0, c_index = 0, r_index = 0; 

    for(int i = 0 ; i < NG ; i++){    
        if(WI[i] != 0) c++;
    }

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] > 0 && WI[i] == 0) {
            D_eig[r_index] = WR[i];

            for(int j = 0 ; j < NG ; j++){
                V_eig[j][r_index] = VR[j][i];
            }

            r_index++;
        }
    }

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] > 0 && WI[i] > 0) {
            D_eig[r_index] = WR[i];
            DI_eig[r_index] = WI[i];

            D_eig[r_index+1] = WR[i];
            DI_eig[r_index+1] = -WI[i];

            for(int j = 0 ; j < NG ; j++){
                V_eig[j][r_index] = VR[j][i];
                VI_eig[j][r_index] = VR[j][i+1];


                V_eig[j][r_index+1] = VR[j][i];
                VI_eig[j][r_index+1] = -VR[j][i+1];
            }

            r_index++;
            r_index++;
        }
    }

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] < 0 && WI[i] == 0) {
            D_eig[r_index] = WR[i];

            for(int j = 0 ; j < NG ; j++){
                V_eig[j][r_index] = VR[j][i];
            }

            r_index++;
        }
    }

    for(int i = 0 ; i < NG ; i++){
        if(WR[i] < 0 && WI[i] > 0) {
            D_eig[r_index] = WR[i];
            DI_eig[r_index] = WI[i];

            D_eig[r_index+1] = WR[i];
            DI_eig[r_index+1] = -WI[i];

            for(int j = 0 ; j < NG ; j++){
                V_eig[j][r_index] = VR[j][i];
                VI_eig[j][r_index] = VR[j][i+1];


                V_eig[j][r_index+1] = VR[j][i];
                VI_eig[j][r_index+1] = -VR[j][i+1];
            }

            r_index++;
            r_index++;
        }
    }

    free(VV);
}
