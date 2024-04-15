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
	 vector<vector<vector<vector<double>>>> &psi_in_y_inv, vector<vector<int>> &zone_config, vector<vector<int>> &zcn){

	int MG = M*G, aux;

	for(int i = 0 ; i  < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M/4 ; k++){
					aux = k + g*M;
					for(int l = 0 ; l < MG ; l++){
						psi_in_x_inv[i][j][aux][l] = V_x[zcn[i][j]][aux][l];
						psi_in_x_inv[i][j][aux + 3*M/4][l] = V_x[zcn[i][j]][aux + 3*M/4][l];
						psi_in_x_inv[i][j][aux + M/4][l] = V_x[zcn[i][j]][aux + M/4][l] * exp(-h_x[i][j]/D_x[zcn[i][j]][l]);
						psi_in_x_inv[i][j][aux + M/2][l] = V_x[zcn[i][j]][aux + M/2][l] * exp(-h_x[i][j]/D_x[zcn[i][j]][l]);
					}
				
					for(int l = 0 ; l < MG ; l++){
						psi_in_y_inv[i][j][aux][l] = V_y[zcn[i][j]][aux][l];
						psi_in_y_inv[i][j][aux + M/4][l] = V_y[zcn[i][j]][aux + M/4][l];
						psi_in_y_inv[i][j][aux + M/2][l] = V_y[zcn[i][j]][aux + M/2][l] * exp(-h_y[i][j]/D_y[zcn[i][j]][l]);
						psi_in_y_inv[i][j][aux + 3*M/4][l] = V_y[zcn[i][j]][aux + 3*M/4][l] * exp(-h_y[i][j]/D_y[zcn[i][j]][l]);
					}
				}
			}

			psi_in_x_inv[i][j] = inverse(psi_in_x_inv[i][j], MG);
			psi_in_y_inv[i][j] = inverse(psi_in_y_inv[i][j], MG);
		}
	}
}

void get_psi_in(const int nxf, const int nyf, const int M, const int G, vector<vector<double>> &D_x,
	 vector<vector<double>> &D_y, vector<vector<vector<double>>> &V_x, vector<vector<vector<double>>> &V_y,
	 vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<vector<double>>>> &psi_in_x, 
	 vector<vector<vector<vector<double>>>> &psi_in_y, vector<vector<int>> &zone_config, vector<vector<int>> &zcn){

	int MG = M*G, aux;

	for(int i = 0 ; i  < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M/4 ; k++){
					aux = k + g*M;
					for(int l = 0 ; l < MG ; l++){
						psi_in_x[i][j][aux][l] = V_x[zcn[i][j]][aux][l] * exp(-h_x[i][j]/D_x[zcn[i][j]][l]);
						psi_in_x[i][j][aux + 3*M/4][l] = V_x[zcn[i][j]][aux + 3*M/4][l] * exp(-h_x[i][j]/D_x[zcn[i][j]][l]);

						psi_in_x[i][j][aux + M/4][l] = V_x[zcn[i][j]][aux + M/4][l];
						psi_in_x[i][j][aux + M/2][l] = V_x[zcn[i][j]][aux + M/2][l];
					}

					for(int l = 0 ; l < MG ; l++){
						psi_in_y[i][j][aux][l] = V_y[zcn[i][j]][aux][l] * exp(-h_y[i][j]/D_y[zcn[i][j]][l]);
						psi_in_y[i][j][aux + M/4][l] = V_y[zcn[i][j]][aux + M/4][l] * exp(-h_y[i][j]/D_y[zcn[i][j]][l]);

						psi_in_y[i][j][aux + M/2][l] = V_y[zcn[i][j]][aux + M/2][l];
						psi_in_y[i][j][aux + 3*M/4][l] = V_y[zcn[i][j]][aux + 3*M/4][l];
					}
				}
			}
		}
	}
}

void get_psi_aux(int nxf, int nyf, vector<vector<vector<vector<double>>>> &psi_x_aux, vector<vector<vector<vector<double>>>> &psi_y_aux,
	vector<vector<vector<vector<double>>>> &psi_in_x, vector<vector<vector<vector<double>>>> &psi_in_y,
	vector<vector<vector<vector<double>>>> &psi_in_x_inv, vector<vector<vector<vector<double>>>> &psi_in_y_inv){
	
	for(int i = 0 ; i  < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			psi_x_aux[i][j] = mult_mat(psi_in_x[i][j], psi_in_x_inv[i][j]);
			psi_y_aux[i][j] = mult_mat(psi_in_y[i][j], psi_in_y_inv[i][j]);
		}
	}
}