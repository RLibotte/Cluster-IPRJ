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

void get_psi_b(int MG, int nxf, int nyf, 
	 vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
	 vector<vector<vector<double>>> &psi_x_b, vector<vector<vector<double>>> &psi_y_b){

	// #pragma omp parallel for collapse(3)
	for(int k = 0 ; k < MG ; k++){
		for(int i = 0 ; i < nyf ; i++){
			for(int j = 0 ; j < nxf ; j++){
				psi_x_b[i][j][k] = (psi_x[i][j][k] + psi_x[i][j+1][k]) * 0.5;
				psi_y_b[i][j][k] = (psi_y[i][j][k] + psi_y[i+1][j][k]) * 0.5;
			}
		}
	}
}