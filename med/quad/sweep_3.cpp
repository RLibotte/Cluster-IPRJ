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
			 vector<vector<vector<double>>> &psi_x_old, 
			 vector<vector<vector<double>>> &psi_y_old, 
			 vector<double> &l_x_temp,
			 vector<double> &l_y_temp,
			 vector<vector<vector<double>>> &q, 
			 vector<double> &b_in_x, 
			 vector<double> &b_in_y, 
			 vector<vector<vector<double>>> &alpha_x, 
			 vector<vector<vector<double>>> &alpha_y,
			 vector<double> &psi_x_p, 
			 vector<double> &psi_y_p, 
			 vector<vector<vector<double>>> &psi_b_x,
			 vector<vector<vector<double>>> &psi_b_y,
			 vector<vector<double>> &h_x, 
			 vector<vector<double>> &h_y, 
			 vector<vector<vector<vector<double>>>> &psi_in_x, 
			 vector<vector<vector<vector<double>>>> &psi_in_y,	
			 vector<vector<vector<vector<double>>>> &psi_in_x_inv, 
			 vector<vector<vector<vector<double>>>> &psi_in_y_inv,
			 vector<vector<vector<double>>> &psi_x_p_out, 
			 vector<vector<vector<double>>> &psi_y_p_out, 
			 vector<double> &psi_x_out, 
			 vector<double> &psi_y_out,
			 vector<vector<vector<double>>> &c1_x, 
			 vector<vector<vector<double>>> &c2_x, 
			 vector<vector<vector<double>>> &c1_y, 
			 vector<vector<vector<double>>> &c2_y, 
			 vector<vector<vector<double>>> &b0_x, 
			 vector<vector<vector<double>>> &b1_x, 
			 vector<vector<vector<double>>> &b2_x, 
			 vector<vector<vector<double>>> &b0_y, 
			 vector<vector<vector<double>>> &b1_y, 
			 vector<vector<vector<double>>> &b2_y,
			 vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z,
			 vector<vector<vector<double>>> &psi_in_x_inv_arr, 
			 vector<vector<vector<double>>> &psi_in_y_inv_arr,
			 vector<vector<vector<double>>> mhx,
			 vector<vector<vector<double>>> nhy){

	int MG = M*G, index;

		for(int i = 0 ; i < nyf ; i++){
			for(int j = nxf - 1 ; j >= 0 ; j--){
			
			quad_l_x(j, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, mhx, nhy, M/4, M/2);
			quad_l_y(j, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, mhx, nhy, M/4, M/2);

			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M/4 ; k++){
					
					index = k + g*M;
					b_in_x[index] = psi_x[i][j][index]     - (b0_x[i][j][index] - b1_x[i][j][index] - b2_x[i][j][index]);
					b_in_y[index] = psi_y[i + 1][j][index] - (b0_y[i][j][index] - b1_y[i][j][index] - b2_y[i][j][index]);

					index += M/4;
					b_in_x[index] = psi_x[i][j + 1][index] - (b0_x[i][j][index] + b1_x[i][j][index] - b2_x[i][j][index]);
					b_in_y[index] = psi_y[i + 1][j][index] - (b0_y[i][j][index] - b1_y[i][j][index] - b2_y[i][j][index]);
			
					index += M/4;
					b_in_x[index] = psi_x[i][j + 1][index] - (b0_x[i][j][index] + b1_x[i][j][index] - b2_x[i][j][index]);
					b_in_y[index] = psi_y[i][j][index]     - (b0_y[i][j][index] + b1_y[i][j][index] - b2_y[i][j][index]);

					index += M/4;
					b_in_x[index] = psi_x[i][j][index]     - (b0_x[i][j][index] - b1_x[i][j][index] - b2_x[i][j][index]);
					b_in_y[index] = psi_y[i][j][index]     - (b0_y[i][j][index] + b1_y[i][j][index] - b2_y[i][j][index]);
				}
			}
			
			alpha_x[i][j] = mult(psi_in_x_inv[i][j], b_in_x);
			alpha_y[i][j] = mult(psi_in_y_inv[i][j], b_in_y);

			psi_x_out = mult(psi_in_x[i][j], alpha_x[i][j]);
			psi_y_out = mult(psi_in_y[i][j], alpha_y[i][j]);

			for(int g = 0 ; g < G ; g++){
				for(int k = M/2 ; k < 3*M/4 ; k++){
					index = k + g*M;

					psi_x[i][j][index] = psi_x_out[index]   + (b0_x[i][j][index] - b1_x[i][j][index] - b2_x[i][j][index]);
					psi_y[i+1][j][index] = psi_y_out[index] + (b0_y[i][j][index] - b1_y[i][j][index] - b2_y[i][j][index]);
				}
			}
		}
	}
}