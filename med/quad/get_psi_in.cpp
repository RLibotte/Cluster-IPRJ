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

void get_psi_in_inv(const int nxf, const int nyf, const int M, const int G, vector<vector<double>> &D_x,
	 vector<vector<double>> &D_y, vector<vector<vector<double>>> &V_x, vector<vector<vector<double>>> &V_y,
	 vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<vector<double>>>> &psi_in_x_inv, 
	 vector<vector<vector<vector<double>>>> &psi_in_y_inv, vector<vector<int>> &zone_config, vector<vector<int>> &zcn,
	 vector<vector<vector<double>>> &psi_in_x_inv_arr, vector<vector<vector<double>>> &psi_in_y_inv_arr){

	int MG = M*G;

	for(int i = 0 ; i  < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M/4 ; k++){
					for(int l = 0 ; l < MG ; l++){
						psi_in_x_inv[i][j][k + M*g][l] = V_x[zcn[i][j]][k + M*g][l];
						psi_in_x_inv[i][j][k + M*g + 3*M/4][l] = V_x[zcn[i][j]][k + M*g + 3*M/4][l];

						psi_in_x_inv[i][j][k + M*g + M/4][l] = V_x[zcn[i][j]][k + M*g + M/4][l] * exp(-h_x[i][j]/D_x[zcn[i][j]][l]);
						psi_in_x_inv[i][j][k + M*g + M/2][l] = V_x[zcn[i][j]][k + M*g + M/2][l] * exp(-h_x[i][j]/D_x[zcn[i][j]][l]);
					}
				
					for(int l = 0 ; l < MG ; l++){
						psi_in_y_inv[i][j][k + M*g][l] = V_y[zcn[i][j]][k + M*g][l];
						psi_in_y_inv[i][j][k + M*g + M/4][l] = V_y[zcn[i][j]][k + M*g + M/4][l];

						psi_in_y_inv[i][j][k + M*g + M/2][l] = V_y[zcn[i][j]][k + M*g + M/2][l] * exp(-h_y[i][j]/D_y[zcn[i][j]][l]);
						psi_in_y_inv[i][j][k + M*g + 3*M/4][l] = V_y[zcn[i][j]][k + M*g + 3*M/4][l] * exp(-h_y[i][j]/D_y[zcn[i][j]][l]);
					}
				}
			}

			psi_in_x_inv[i][j] = inverse(psi_in_x_inv[i][j], MG);
			psi_in_y_inv[i][j] = inverse(psi_in_y_inv[i][j], MG);

			for(int k = 0 ; k < MG ; k++){
				for(int l = 0 ; l < MG ; l++){
					psi_in_x_inv_arr[i][j][k * MG + l] = psi_in_x_inv[i][j][k][l];
					psi_in_y_inv_arr[i][j][k * MG + l] = psi_in_y_inv[i][j][k][l];
				}
			}
		}
	}
}

void get_psi_in(const int nxf, const int nyf, const int M, const int G, vector<vector<double>> &D_x,
	 vector<vector<double>> &D_y, vector<vector<vector<double>>> &V_x, vector<vector<vector<double>>> &V_y,
	 vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<vector<double>>>> &psi_in_x, 
	 vector<vector<vector<vector<double>>>> &psi_in_y, vector<vector<int>> &zone_config, vector<vector<int>> &zcn){

	int MG = M*G;

	for(int i = 0 ; i  < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M/4 ; k++){
					for(int l = 0 ; l < MG ; l++){
						psi_in_x[i][j][k + g*M][l] = V_x[zcn[i][j]][k + g*M][l] * exp(-h_x[i][j]/D_x[zcn[i][j]][l]);
						psi_in_x[i][j][k + g*M + 3*M/4][l] = V_x[zcn[i][j]][k + g*M + 3*M/4][l] * exp(-h_x[i][j]/D_x[zcn[i][j]][l]);

						psi_in_x[i][j][k + g*M + M/4][l] = V_x[zcn[i][j]][k + g*M + M/4][l];
						psi_in_x[i][j][k + g*M + M/2][l] = V_x[zcn[i][j]][k + g*M + M/2][l];
					}

					for(int l = 0 ; l < MG ; l++){
						psi_in_y[i][j][k + g*M][l] = V_y[zcn[i][j]][k + g*M][l] * exp(-h_y[i][j]/D_y[zcn[i][j]][l]);
						psi_in_y[i][j][k + g*M + M/4][l] = V_y[zcn[i][j]][k + g*M + M/4][l] * exp(-h_y[i][j]/D_y[zcn[i][j]][l]);

						psi_in_y[i][j][k + g*M + M/2][l] = V_y[zcn[i][j]][k + g*M + M/2][l];
						psi_in_y[i][j][k + g*M + 3*M/4][l] = V_y[zcn[i][j]][k + g*M + 3*M/4][l];
					}
				}
			}
		}
	}
}