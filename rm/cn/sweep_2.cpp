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

void sweep_2(const int &M, 
			 const int &G, 
			 int &nxf, 
			 int &nyf, 
			 vector<vector<vector<double>>> &l_x,
	 		 vector<vector<vector<double>>> &l_y, 
			 vector<vector<vector<double>>> &psi_x, 
			 vector<vector<vector<double>>> &psi_y, 
			 vector<vector<vector<double>>> &mi_h, 
			 vector<vector<vector<double>>> &n_h, 
			 vector<double> &l_x_temp,
			 vector<double> &l_y_temp,
			 vector<vector<vector<double>>> &q, 
			 vector<double> &b_in_x, 
			 vector<double> &b_in_y, 
			 vector<vector<double>> &h_x, 
			 vector<vector<double>> &h_y, 
			 vector<vector<vector<vector<double>>>> &psi_iter_x, 
			 vector<vector<vector<vector<double>>>> &psi_iter_y, 
			 vector<vector<vector<vector<double>>>> &psi_aux_x, 
			 vector<vector<vector<vector<double>>>> &psi_aux_y){

	int MG = M*G, index;

	for(int i = nyf - 1 ; i >= 0 ; i--){
		for(int j = nxf - 1 ; j >= 0 ; j--){
			for(int g = 0 ; g < G ; ++g){
				for(int k = 0 ; k < M/4 ; k++){
					index = k + g*M;

					l_x[i][j][index] = mi_h[i][j][k] * (psi_x[i][j+1][index] - psi_x[i][j][index]);
					l_y[i][j][index] =  n_h[i][j][k] * (psi_y[i][j][index] - psi_y[i+1][j][index]);
				}
			}

			for(int g = 0 ; g < G ; ++g){
				for(int k = 0 ; k < M ; k++){
					index = k + g*M;

					l_x_temp[index] = -q[i][j][g] + l_x[i][j][index];
					l_y_temp[index] = -q[i][j][g] + l_y[i][j][index];
				}
			}
			
			get_b_in(i, j, G, M, b_in_x, b_in_y, psi_x, psi_y);

			mult_part(psi_iter_x[i][j], b_in_x, b_in_x, M/4, M/2, G, M);
			mult_part(psi_aux_x[i][j], l_y_temp, l_y_temp, M/4, M/2, G, M);

			mult_part(psi_iter_y[i][j], b_in_y, b_in_y, M/4, M/2, G, M);
			mult_part(psi_aux_y[i][j], l_x_temp, l_x_temp, M/4, M/2, G, M);

			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M/4 ; k++){
					index =  k + g*M + M/4;

					psi_x[i][j][index] = b_in_x[index] + l_y_temp[index];
					psi_y[i][j][index] = b_in_y[index] + l_x_temp[index];
				}
			}
		}
	}

}