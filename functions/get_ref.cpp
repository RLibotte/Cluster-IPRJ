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

void get_ref_l(const int M, const int G, int nxf, int nyf, 
	vector<vector<vector<double>>> &psi_x){
	for(int g = 0 ; g < G ; g++){
		for(int k = M/4 ; k < M/2 ; k++){
			for(int i = 0 ; i < nyf ; i++){
				psi_x[i][nxf][k + M*g] = psi_x[i][nxf][k - M/4 + M*g];
			}
		}
	}

	for(int g = 0 ; g < G ; g++){
		for(int k = M/2 ; k < 3*M/4 ; k++){
			for(int i = 0 ; i < nyf ; i++){
				psi_x[i][nxf][k + M*g] = psi_x[i][nxf][k + M/4 + M*g];
			}
		}
	}
}

void get_ref_o(const int M, const int G, int nxf, int nyf, 
	vector<vector<vector<double>>> &psi_x){
	for(int g = 0 ; g < G ; g++){
		for(int k = 0 ; k < M/4 ; k++){
			for(int i = 0 ; i < nyf ; i++){
				psi_x[i][0][k + g*M] = psi_x[i][0][k + M/4 + g*M];
			}
		}
	}
	for(int g = 0 ; g < G ; g++){
		for(int k = 3*M/4 ; k < M ; k++){
			for(int i = 0 ; i < nyf ; i++){
				psi_x[i][0][k + g*M] = psi_x[i][0][k - M/4 + g*M];
			}
		}
	}
}

void get_ref_n(const int M, const int G, int nxf, int nyf, 
	vector<vector<vector<double>>> &psi_y){
	for(int g = 0 ; g < G ; g++){
		for(int k = M/2 ; k < 3*M/4 ; k++){
			for(int i = 0 ; i < nxf ; i++){
				psi_y[0][i][k + M*g] = psi_y[0][i][k - M/4 + M*g];
			}
		}
	}

	for(int g = 0 ; g < G ; g++){
		for(int k = 3*M/4 ; k < M ; k++){
			for(int i = 0 ; i < nxf ; i++){
				psi_y[0][i][k + M*g] = psi_y[0][i][k - 3*M/4 + M*g];
			}
		}
	}
}

void get_ref_s(const int M, const int G, int nxf, int nyf, 
	vector<vector<vector<double>>> &psi_y){
	for(int g = 0 ; g < G ; g++){
		for(int k = 0 ; k < M/4 ; k++){
			for(int i = 0 ; i < nxf ; i++){
				psi_y[nyf][i][k + g*M] = psi_y[nyf][i][k + 3*M/4 + g*M];
			}
		}
	}
	for(int g = 0 ; g < G ; g++){
		for(int k = M/4 ; k < M/2 ; k++){
			for(int i = 0 ; i < nxf ; i++){
				psi_y[nyf][i][k + g*M] = psi_y[nyf][i][k + M/4 + g*M];
			}
		}
	}
}