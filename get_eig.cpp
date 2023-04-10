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

void get_eig(int R, int Z, int G, int N, std::vector<std::vector<double>> &sigmat, 
	 std::vector<std::vector<std::vector <double>>> &sigmas0, std::vector<std::vector<std::vector <double>>> &sigmas1, 
     std::vector<double> &mi, std::vector<double> &w, std::vector<int> &zone_config, 
     std::vector<std::vector<double>> &ni, std::vector<std::vector<double>> &ni_i,
     std::vector<std::vector<std::vector <double>>> &V, std::vector<std::vector<std::vector <double>>> &VI){

	const int NG = N*G;
	int r_index = 0, c_index = 0;

    D.resize(NG);

	A.resize(NG, std::vector<double> (NG));
	V_temp.resize(NG, std::vector<double> (NG));

	A_p.resize(NG, std::vector<std::vector<double>> (NG, std::vector<double> (R)));

	for (int i = 0; i < Z; i++) {
        for (int j = 0; j < G; j++) {
            for (int k = 0; k < N; k++) {
                for (int l = 0; l < G; l++) {
                    for (int m = 0; m < N; m++) {
                        if (r_index == c_index) {
                            A[r_index][c_index] = (sigmat[i][j] / mi[k]) - (sigmas0[i][l][j] * w[m]) / (2.0 * mi[k]) - (1.5 * sigmas1[i][l][j] * mi[m] * w[m]);
                            A_p[r_index][c_index][i] = (sigmat[i][j]) - (sigmas0[i][l][j] * w[m]) / (2.0) - (1.5 * sigmas1[i][l][j] * mi[m] * w[m]);
                            c_index++;
                        } else {
                            A[r_index][c_index] = -(sigmas0[i][l][j] * w[m]) / (2.0 * mi[k]) - (1.5 * sigmas1[i][l][j] * mi[m] * w[m]);
                            A_p[r_index][c_index][i] = -(sigmas0[i][l][j] * w[m]) / (2.0) - (1.5 * sigmas1[i][l][j] * mi[m] * w[m]);
                            c_index++;
                        }
                    }
                }
                r_index++;
                c_index = 0;
            }
        }
        
        eigen(A, NG, D_eig, DI_eig, V_eig, VI_eig);

        for(int l = 0 ; l < R ; l++){
            if(zone_config[l] == (i+1)){
                for(int j = 0 ; j < NG ; j++){
                    for(int k = 0 ; k < NG ; k++){
                        V[l][j][k] = V_eig[j][k];
                        VI[l][j][k] = VI_eig[j][k];
                    }
                    ni[j][l] = D_eig[j];
                    ni_i[j][l] = DI_eig[j];
                }   
            }
        }
        r_index = 0;
    }
}