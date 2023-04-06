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

tuple <vector<vector<double>>, vector<vector<vector<double>>>, vector<vector<vector<double>>>> 
	get_eig_h(int R, int Z, int G, int N, vector<vector<double>> sigmat, 
	vector<vector<vector <double>>> sigmas0, vector<vector<vector <double>>> sigmas1, 
    vector<double> mi, vector<double> w, vector<int> zone_config){

	const int NG = N*G;
	int r_index = 0, c_index = 0;

    vector<double> D (NG);

	vector<vector<double>> A (NG/2, vector<double> (NG/2));
	vector<vector<double>> ni (NG, vector<double> (R));
	vector<vector<double>> V_temp (NG, vector<double> (NG));

	vector<vector<vector<double>>> V(R, vector<vector<double>> (NG, vector<double> (NG)));
	vector<vector<vector<double>>> A_p(NG, vector<vector<double>> (NG, vector<double> (R)));

	for (int i = 0; i < Z; i++) {
        for (int j = 0; j < G; j++) {
            for (int k = 0; k < N; k++) {
                for (int l = 0; l < G; l++) {
                    for (int m = 0; m < N; m++) {
                        if (r_index == c_index) {
                            A_p[r_index][c_index][i] = (sigmat[i][j]) - (sigmas0[i][l][j] * w[m]) / (2.0);// - (1.5 * sigmas1[i][l][j] * mi[m] * w[m]);
                            c_index++;
                        } else {
                            A_p[r_index][c_index][i] = -(sigmas0[i][l][j] * w[m]) / (2.0);// - (1.5 * sigmas1[i][l][j] * mi[m] * w[m]);
                            c_index++;
                        }
                    }
                }
                r_index++;
                c_index = 0;
            }
        }
        r_index = 0;
    }

	for (int i = 0; i < Z; i++) {
        for (int j = 0; j < G; j++) {
            for (int k = 0; k < N/2; k++) {
                for (int l = 0; l < G; l++) {
                    for (int m = 0; m < N/2; m++) {
                        if (r_index == c_index) {
                            A[r_index][c_index] = ((sigmat[i][j] / mi[k]) * (sigmat[i][j] / mi[k])) - ((sigmat[i][j] * sigmas0[i][l][j] * w[m]) / (mi[k] * mi[k]));// - (1.5 * sigmas1[i][l][j] * mi[m] * w[m]);
                            c_index++;
                        } else {
                            A[r_index][c_index] = -(sigmat[i][j] * sigmas0[i][l][j] * w[m]) / (mi[k] * mi[k]);// - (1.5 * sigmas1[i][l][j] * mi[m] * w[m]);
                            c_index++;
                        }
                    }
                }
                r_index++;
                c_index = 0;
            }
        }

        auto [V_temp, D] = eigen(A, NG/2);

        for(int l = 0 ; l < R ; l++){
            if(zone_config[l] == (i+1)){
                for(int g = 0 ; g < G ; g++){
                    for(int j = 0 ; j < N/2 ; j++){
                        ni[j + N*g][l] = 1.0/sqrt(abs(D[j + N*g/2]));
                        ni[j+N/2 + N*g][l] = -1.0/sqrt(abs(D[j + N*g/2]));
                    }
                }
                for(int g = 0 ; g < G ; g++){
                    for(int j = 0 ; j < N/2 ; j++){
                        for(int h = 0 ; h < G ; h++){
                            for(int k = 0 ; k < N/2 ; k++){
                                V[l][j + N*g][k + N*h] = (V_temp[j + N*g/2][k + N*h/2]/2.0)*(1.0 + mi[j]/(ni[k + N*h][l]*sigmat[i][g]));
                                V[l][j + N/2 + N*g][k + N*h] = (V_temp[j + N*g/2][k + N*h/2]/2.0)*(1.0 - mi[j]/(ni[k + N*h][l]*sigmat[i][g]));

                                V[l][j + N*g][k + N/2 + N*h] = (V_temp[j + N*g/2][k + N*h/2]/2.0)*(1.0 - mi[j]/(ni[k + N*h][l]*sigmat[i][g]));
                                V[l][j + N/2 + N*g][k + N/2 + N*h] = (V_temp[j + N*g/2][k + N*h/2]/2.0)*(1.0 + mi[j]/(ni[k + N*h][l]*sigmat[i][g]));
                            }
                        }   
                    }
                }
            }
        }
         
        r_index = 0;
    }

    return {ni, V, A_p};
}