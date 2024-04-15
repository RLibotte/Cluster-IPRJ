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

	 get_matrix(const int &M, const int &G, int &nxf, int &nyf, int x_size, int y_size,
	 			vector<vector<double>> &D_x, vector<vector<double>> &D_y,
	 			vector<vector<vector<double>>> &V_x, vector<vector<vector<double>>> &V_y,
	 			vector<vector<double>> &h_x, vector<vector<double>> &h_y,vector<vector<double>> &h_x_n, vector<vector<double>> &h_y_n,
	 			vector<vector<int>> &nodes_x, vector<vector<int>> &nodes_y, vector<vector<vector<double>>> &A_p_inv,
	 			vector<vector<int>> &zone_config, vector<vector<int>> &zcn, vector<double> &mi, vector<double> &n){

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
							temp_in_y(MG , vector<double> (MG));

	vector<vector<vector<double>>>  psi_x (nyf, vector<vector<double>> (nxf + 1, vector <double>(MG,0))),
									psi_x_old (nyf, vector<vector<double>> (nxf + 1, vector <double>(MG,0))), 
								    psi_y (nyf + 1, vector<vector<double>> (nxf, vector <double>(MG,0))),
								    psi_b_x (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))), 
								    psi_b_y (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))),
								    alpha_x (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))), 
								    alpha_y (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))),
								    l_x (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))), 
								    l_y (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))),
								    l_x_old (nyf, vector<vector<double>> (nxf, vector <double>(MG,0))),
									scalar_flux_x (nyf, vector<vector<double>> (nxf, vector <double>(G,0))),
									scalar_flux_x_out (y_size, vector<vector<double>> (x_size, vector <double>(G,0))),
									psi_x_p_out (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
									psi_y_p_out (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								    mhx (nyf, vector<vector<double>> (nxf, vector <double>(M,0))),
								    nhy (nyf, vector<vector<double>> (nxf, vector <double>(M,0)));


	vector<vector<vector<vector<double>>>> A_p_lu (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_x_p (y_size, vector<vector<vector<double>>> (x_size, vector<vector<double>>(MG, vector<double> (MG)))),
									 	   psi_in_y_p (y_size, vector<vector<vector<double>>> (x_size, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_x_lu (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_y_lu (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_x (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_y (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_iter_x (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_iter_y (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_x_inv (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_in_y_inv (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_x_aux (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG)))),
										   psi_y_aux (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(MG, vector<double> (MG))));


	for(int i = 0 ; i < nyf ; ++i){
		for(int j = 0 ; j < nxf ; ++j){
			for(int k = 0 ; k < M ; k++){
				mhx[i][j][k] = mi[k]/h_x[i][j];
				nhy[i][j][k] = n[k]/h_y[i][j];
			}
		}
	}

	index_x = 0, index_y = 0;
	for(int l = 0 ; l < y_size ; l++){
		for(int k = 0 ; k < x_size ; k++){
			for(int a = 0 ; a < nodes_y[l][k] ; ++a){
				for(int b = 0 ; b < nodes_x[l][k] ; ++b){
					A_p_lu[a + index_y][b + index_x] = A_p_inv[zone_config[l][k]];
				}
			}
			index_x += nodes_x[l][k];
		}
		index_x = 0;
		index_y += nodes_y[l][0];
	}

	return {b_in_x, b_in_y, psi_x_p, psi_y_p, l_x_temp, l_y_temp, psi_x_out, psi_y_out,
			psi_x, psi_x_old, psi_y, psi_b_x, psi_b_y, alpha_x, alpha_y, l_x, l_y, l_x_old, scalar_flux_x, scalar_flux_x_out, psi_x_p_out,
			psi_y_p_out, mhx, nhy,
			A_p_lu, psi_in_x, psi_in_y, psi_in_x_inv, psi_in_y_inv, psi_iter_x, psi_iter_y, psi_x_aux, psi_y_aux};
}