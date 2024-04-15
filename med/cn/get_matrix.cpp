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


tuple<vector<double>,
		vector<double>,
		vector<double>,
		vector<double>,
		vector<double>,
		vector<double>,
		vector<double>,
		vector<double>,
		vector<vector<double>>,
		vector<vector<double>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<double>>>,
		vector<vector<vector<vector<double>>>>,
		vector<vector<vector<vector<double>>>>,
		vector<vector<vector<vector<double>>>>,
		vector<vector<vector<vector<double>>>>,
		vector<vector<vector<vector<double>>>>,
		vector<vector<vector<vector<double>>>>,
		vector<vector<vector<vector<double>>>>,
		vector<vector<vector<vector<double>>>>,
		vector<vector<vector<vector<double>>>>>

	 get_matrix(const int M, const int G, int nxf, int nyf, int x_size, int y_size,
	 			vector<vector<double>> D_x, vector<vector<double>> D_y,
	 			vector<vector<vector<double>>> V_x, vector<vector<vector<double>>> V_y,
	 			vector<vector<double>> h_x, vector<vector<double>> h_y,vector<vector<double>> h_x_n, vector<vector<double>> h_y_n,
	 			vector<vector<int>> nodes_x, vector<vector<int>> nodes_y, vector<vector<vector<double>>> A_p_inv,
	 			vector<vector<int>> zone_config, vector<vector<int>> zcn){

	const int MG = M * G;

	int index_x = 0, index_y = 0, r_index = 0, c_index = 0, z = 0;

	double dev_x = 0;

	vector<double> b_in_x(MG), 
				   b_in_y(MG), 
				   psi_x_p (MG), 
				   psi_y_p (MG), 
				   l_x_temp (MG), 
				   l_y_temp (MG),
				   psi_x_out (MG),
				   psi_y_out (MG);

	vector<vector <double>> temp_x(MG , vector<double> (MG)),
							temp_in_x(MG , vector<double> (MG)),
							temp_in_y(MG , vector<double> (MG)),
							SS_x(G , vector<double> (2)),
							SS_y(G , vector<double> (2));

	vector<vector<vector<double>>>  psi_x (nyf, vector<vector<double>> (nxf + 1, vector <double>(MG,0))), 
								    psi_y (nyf + 1, vector<vector<double>> (nxf, vector <double>(MG,0))),
								    psi_b_x (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))), 
								    psi_b_y (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))),
								    alpha_x (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))), 
								    alpha_y (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))),
								    l_x (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))), 
								    l_y (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))),
									scalar_flux_x (nyf, vector<vector<double>> (nxf, vector <double>(G,0))),
 									scalar_flux_x_old (nyf, vector<vector<double>> (nxf, vector <double>(G,0))),
									scalar_flux_x_out (y_size, vector<vector<double>> (x_size, vector <double>(G,0))),
									psi_x_p_out (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
									psi_y_p_out (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
									psi_x_old (nyf, vector<vector<double>> (nxf + 1, vector <double>(MG))),
									psi_y_old (nyf + 1, vector<vector<double>> (nxf, vector <double>(MG)));


	vector<vector<vector<vector<double>>>> A_p_lu (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_x_p (y_size, vector<vector<vector<double>>> (x_size, vector<vector<double>>(MG, vector<double> (MG)))),
									 	   psi_in_y_p (y_size, vector<vector<vector<double>>> (x_size, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_x_lu (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_y_lu (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_b_x_aux (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_b_y_aux (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_x (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_y (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_x_inv (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_y_inv (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG))));


	for(int i = 0 ; i < nyf ; ++i){
		for(int j = 0 ; j < nxf ; ++j){
			for(int k = 0 ; k < MG ; ++k){
				for(int l = 0 ; l < MG ; ++l){
					if(D_x[zcn[i][j]][l] > 0){
						psi_b_x_aux[i][j][k][l] = D_x[zcn[i][j]][l] * V_x[zcn[i][j]][k][l] * (exp(-h_x[i][j]/D_x[zcn[i][j]][l]) - 1.0)/h_x[i][j];
					} else {
						psi_b_x_aux[i][j][k][l] = D_x[zcn[i][j]][l] * V_x[zcn[i][j]][k][l] * (1.0 - exp(h_x[i][j]/D_x[zcn[i][j]][l]))/h_x[i][j];
					}
					if(D_y[zcn[i][j]][l] > 0){
						psi_b_y_aux[i][j][k][l] = D_y[zcn[i][j]][l] * V_y[zcn[i][j]][k][l] * (exp(-h_y[i][j]/D_y[zcn[i][j]][l]) - 1.0)/h_y[i][j];
					} else {
						psi_b_y_aux[i][j][k][l] = D_y[zcn[i][j]][l] * V_y[zcn[i][j]][k][l] * (1.0 - exp(h_y[i][j]/D_y[zcn[i][j]][l]))/h_y[i][j];	
					}
				}
			}
		}
	}

	index_x = 0, index_y = 0;
	for(int l = 0 ; l < y_size ; l++){
		for(int k = 0 ; k < x_size ; k++){

			psi_in_x_p[l][k] = V_x[zone_config[l][k]];
			psi_in_y_p[l][k] = V_y[zone_config[l][k]];

			for(int i = 0 ; i < M/4 ; ++i){
				for(int j = 0 ; j < MG ; ++j){
					for(int g = 0 ; g < G ; ++g){
						if(D_x[zone_config[l][k]][j] > 0){
							psi_in_x_p[l][k][i + M/4 + g*M][j] = psi_in_x_p[l][k][i + M/4 + g*M][j] * exp(-h_x_n[l][k]/(D_x[zone_config[l][k]][j] * nodes_x[l][k]));
							psi_in_x_p[l][k][i + M/2 + g*M][j] = psi_in_x_p[l][k][i + M/2 + g*M][j] * exp(-h_x_n[l][k]/(D_x[zone_config[l][k]][j] * nodes_x[l][k]));
						} else {
							psi_in_x_p[l][k][i + g*M][j] = psi_in_x_p[l][k][i + g*M][j] * exp(h_x_n[l][k]/(D_x[zone_config[l][k]][j] * nodes_x[l][k]));
							psi_in_x_p[l][k][i + 3*M/4 + g*M][j] = psi_in_x_p[l][k][i + 3*M/4 + g*M][j] * exp(h_x_n[l][k]/(D_x[zone_config[l][k]][j] * nodes_x[l][k]));
						}

						if(D_y[zone_config[l][k]][j] > 0){
							psi_in_y_p[l][k][i + M/2 + g*M][j] = psi_in_y_p[l][k][i + M/2 + g*M][j] * exp(-h_y_n[l][k]/(D_y[zone_config[l][k]][j] * nodes_y[l][k]));
							psi_in_y_p[l][k][i + 3*M/4 + g*M][j] = psi_in_y_p[l][k][i + 3*M/4 + g*M][j] * exp(-h_y_n[l][k]/(D_y[zone_config[l][k]][j] * nodes_y[l][k]));
						} else {
							psi_in_y_p[l][k][i + g*M][j] = psi_in_y_p[l][k][i + g*M][j] * exp(h_y_n[l][k]/(D_y[zone_config[l][k]][j] * nodes_y[l][k]));
							psi_in_y_p[l][k][i + M/4 + g*M][j] = psi_in_y_p[l][k][i + M/4 + g*M][j] * exp(h_y_n[l][k]/(D_y[zone_config[l][k]][j] * nodes_y[l][k]));
							
						}
					}
				}
			}

			temp_x = A_p_inv[zone_config[l][k]], MG;
			temp_in_x = inverse(psi_in_x_p[l][k], MG);
			temp_in_y = inverse(psi_in_y_p[l][k], MG);

			for(int a = 0 ; a < nodes_y[l][k] ; ++a){
				for(int b = 0 ; b < nodes_x[l][k] ; ++b){
					A_p_lu[a + index_y][b + index_x] = temp_x;
					psi_in_x_lu[a + index_y][b + index_x] = temp_in_x;
					psi_in_y_lu[a + index_y][b + index_x] = temp_in_y;
				}
			}
			index_x += nodes_x[l][k];
		}
		index_x = 0;
		index_y += nodes_y[l][0];
	}

	return {b_in_x, b_in_y, psi_x_p, psi_y_p, l_x_temp, l_y_temp, psi_x_out, psi_y_out,
			SS_x, SS_y, 
			psi_x, psi_y, psi_b_x, psi_b_y, alpha_x, alpha_y, l_x, l_y, scalar_flux_x, scalar_flux_x_old, scalar_flux_x_out, psi_x_p_out,
			psi_y_p_out, psi_x_old, psi_y_old,
			A_p_lu, psi_in_x_lu, psi_in_y_lu, psi_b_x_aux, psi_b_y_aux, psi_in_x, psi_in_y, psi_in_x_inv, psi_in_y_inv};
}