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

void get_psi_b(const int nxf, const int nyf, const int M, const int G, vector<vector<vector<double>>> &psi_b_x,
			   vector<vector<vector<double>>> &psi_b_y, vector<vector<double>> D_x, vector<vector<double>> D_y,
			   vector<vector<vector<double>>> alpha_x, vector<vector<vector<double>>> alpha_y, vector<vector<vector<double>>> V_x,
			   vector<vector<vector<double>>> V_y, vector<vector<double>> h_x, vector<vector<double>> h_y, 
			   vector<vector<vector<double>>> psi_x_p_out, vector<vector<vector<double>>> psi_y_p_out, vector<vector<int>> zcn){


	for(int i = 0 ; i < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M ; k++){
					psi_b_x[i][j][k + g*M] = 0;
					psi_b_y[i][j][k + g*M] = 0;
					for(int l = 0 ; l < M*G ; l++){
							psi_b_x[i][j][k + g*M] = psi_b_x[i][j][k + g*M] + D_x[zcn[i][j]][l] * alpha_x[i][j][l] * V_x[zcn[i][j]][k + g*M][l] * (exp(-h_x[i][j]/D_x[zcn[i][j]][l]) - 1);
							psi_b_y[i][j][k + g*M] = psi_b_y[i][j][k + g*M] + D_y[zcn[i][j]][l] * alpha_y[i][j][l] * V_y[zcn[i][j]][k + g*M][l] * (exp(-h_y[i][j]/D_y[zcn[i][j]][l]) - 1);
					}
					psi_b_x[i][j][k + g*M] = - psi_b_x[i][j][k + g*M] / h_x[i][j] + psi_x_p_out[i][j][k + g*M];
					psi_b_y[i][j][k + g*M] = - psi_b_y[i][j][k + g*M] / h_y[i][j] + psi_y_p_out[i][j][k + g*M];
				}
			}
		}
	}
}