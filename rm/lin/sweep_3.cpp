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

void sweep_3(const int &M, 
			 const int &G, 
			 int &nxf, 
			 int &nyf, 
			 vector<vector<vector<double>>> &l_x,
	 		 vector<vector<vector<double>>> &l_y, 
	 		 vector<double> &mi, 
	 		 vector<double> &n, 
	 		 vector<double> &w,
			 vector<vector<vector<double>>> &psi_x, 
			 vector<vector<vector<double>>> &psi_y, 
			 vector<vector<vector<double>>> &q,
			 vector<double> &b_in_x, 
			 vector<double> &b_in_y,
			 vector<vector<double>> &h_x, 
			 vector<vector<double>> &h_y,
			 vector<double> &psi_x_out, 
			 vector<double> &psi_y_out,
			 vector<vector<vector<double>>> &c1_x, 
			 vector<vector<vector<double>>> &c1_y, 
			 vector<vector<vector<double>>> &b0_x, 
			 vector<vector<vector<double>>> &b1_x, 
			 vector<vector<vector<double>>> &b0_y, 
			 vector<vector<vector<double>>> &b1_y,
			 vector<vector<int>> &zcn, 
			 vector<vector<vector<double>>> &A_p_inv_z,
			 vector<vector<vector<double>>> &mhx, 
			 vector<vector<vector<double>>> &nhy,
			 vector<vector<vector<vector<double>>>> &psi_x_aux,
			 vector<vector<vector<vector<double>>>> &psi_y_aux){

		int index;

		for(int i = 0 ; i < nyf ; i++){
			for(int j = nxf - 1 ; j >= 0 ; j--){
			
			quad_l_x(j, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, b0_x, b1_x, psi_x, psi_y, M/4, M/2, mhx, nhy);
			quad_l_y(j, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, b0_y, b1_y, psi_x, psi_y, M/4, M/2, mhx, nhy);

			get_b_in(i, j, G, M, b_in_x, b_in_y, b0_x, b1_x, b0_y, b1_y, psi_x, psi_y);

			mult_part(psi_x_aux[i][j], b_in_x, psi_x_out, M/2, 3*M/4, G, M);
			mult_part(psi_y_aux[i][j], b_in_y, psi_y_out, M/2, 3*M/4, G, M);

			for(int g = 0 ; g < G ; g++){
				for(int k = M/2 ; k < 3*M/4 ; k++){
					index = k + g*M;

					psi_x[i][j][index] = psi_x_out[index]   + (b0_x[i][j][index] - b1_x[i][j][index]);
					psi_y[i+1][j][index] = psi_y_out[index] + (b0_y[i][j][index] - b1_y[i][j][index]);
				}
			}
		}
	}
}