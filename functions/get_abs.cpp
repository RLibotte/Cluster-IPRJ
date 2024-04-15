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

tuple<vector<vector<vector<double>>>, vector<vector<vector<double>>>>
			get_abs(const int G, int x_size, int y_size, int nxf, int nyf,
			vector<vector<double>> sigmat_z, vector<vector<vector<double>>> sigmas0_z,
			vector<vector<double>> h_x_n, vector<vector<double>> h_y_n,
			vector<vector<double>> h_x, vector<vector<double>> h_y,
			vector<vector<vector<double>>> scalar_flux_x,vector<vector<vector<double>>> scalar_flux_x_out,
			vector<vector<int>> zone_config, vector<vector<int>> zcn){

	vector<vector<vector<double>>>  t_s0 (nyf, vector<vector<double>> (nxf, vector <double>(G,0))),
								    t_abs_x (nyf, vector<vector<double>> (nxf, vector <double>(G,0))),
									t_s0_n (y_size, vector<vector<double>> (x_size, vector <double>(G,0))),
									t_abs_x_out (y_size, vector<vector<double>> (x_size, vector <double>(G,0)));

	for(int i = 0 ; i < y_size ; i++){
    	for(int j = 0 ; j < x_size ; j++){
    		for(int g = 0 ; g < G ; g++){
    			for(int k = 0 ; k < G ; k++){
    				t_s0_n[i][j][g] = t_s0_n[i][j][g] + sigmas0_z[zone_config[i][j]][g][k];
				}
    		}
    	}
    }

	for(int i = 0 ; i < y_size ; i++){
		for(int j = 0 ; j < x_size ; j++){
			for(int g = 0 ; g < G ; g++){
				t_abs_x_out[i][j][g] = scalar_flux_x_out[i][j][g] * (sigmat_z[zone_config[i][j]][g] - t_s0_n[i][j][g]) * (h_x_n[i][j] * h_y_n[i][j]);
			}
		}
	}

	for(int i = 0 ; i < nyf ; i++){
    	for(int j = 0 ; j < nxf ; j++){
    		for(int g = 0 ; g < G ; g++){
    			for(int k = 0 ; k < G ; k++){
    				t_s0[i][j][g] = t_s0[i][j][g] + sigmas0_z[zcn[i][j]][k][g];
    			}
    		}
    	}
    }

	for(int i = 0 ; i < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; g++){
				t_abs_x[i][j][g] = scalar_flux_x[i][j][g] * (sigmat_z[zcn[i][j]][g] - t_s0[i][j][g]) * (h_x[i][j] * h_y[i][j]);
			}
		}
	}

	return {t_abs_x_out, t_abs_x};
}