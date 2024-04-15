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

void sweep_4(const int &M, 
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
			 vector<vector<vector<vector<double>>>> &A_p_lu,
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
			 vector<vector<vector<vector<double>>>> &psi_in_x_lu, 
			 vector<vector<vector<vector<double>>>> &psi_in_y_lu, 
			 vector<vector<vector<vector<double>>>> &psi_in_x, 
			 vector<vector<vector<vector<double>>>> &psi_in_y,	
			 vector<vector<vector<vector<double>>>> &psi_in_x_inv, 
			 vector<vector<vector<vector<double>>>> &psi_in_y_inv,
			 vector<vector<double>> &SS_x, 
			 vector<vector<double>> &SS_y,
			 vector<vector<vector<double>>> &psi_x_p_out, 
			 vector<vector<vector<double>>> &psi_y_p_out, 
			 vector<double> &psi_x_out, 
			 vector<double> &psi_y_out,
			 vector<vector<vector<vector<double>>>> &psi_b_x_aux,
			 vector<vector<vector<vector<double>>>> &psi_b_y_aux){

	int MG = M*G;

		for(int i = 0 ; i < nyf ; i++){
			for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; ++g){
				for(int k = M/2 ; k < 3*M/4 ; ++k){
					l_x[i][j][k + g*M] = (mi[k]/h_x[i][j]) * (psi_x[i][j+1][k + g*M] - psi_x[i][j][k + g*M]);
					l_y[i][j][k + g*M] =  (n[k]/h_y[i][j]) * (psi_y[i][j][k + g*M] - psi_y[i+1][j][k + g*M]);
				}
			}

			for(int g = 0 ; g < G ; ++g){
				for(int k = 0 ; k < M ; k++){
					l_x_temp[k + g*M] = q[i][j][g] - l_x[i][j][k + g*M];
					l_y_temp[k + g*M] = q[i][j][g] - l_y[i][j][k + g*M];
				}
			}

			psi_x_p = mult(A_p_lu[i][j], l_y_temp);
			psi_y_p = mult(A_p_lu[i][j], l_x_temp);

			#pragma omp parallel for collapse(2)
			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M/4 ; k++){
					b_in_x[k + M*g] = psi_x[i][j][k + M*g] - psi_x_p[k + M*g];
					b_in_x[k + M*g + M/4] = (psi_x[i][j + 1][k + M*g + M/4] - psi_x_p[k + M*g + M/4]);
					b_in_x[k + M*g + M/2] = (psi_x[i][j + 1][k + M*g + M/2] - psi_x_p[k + M*g + M/2]);
					b_in_x[k + M*g + 3*M/4] = (psi_x[i][j][k + M*g + 3*M/4] - psi_x_p[k + M*g + 3*M/4]);

					b_in_y[k + M*g] = (psi_y[i + 1][j][k + M*g] - psi_y_p[k + M*g]);
					b_in_y[k + M*g + M/4] = (psi_y[i + 1][j][k + M*g + M/4] - psi_y_p[k + M*g + M/4]);
					b_in_y[k + M*g + M/2] = (psi_y[i][j][k + M*g + M/2] - psi_y_p[k + M*g + M/2]);
					b_in_y[k + M*g + 3*M/4] = (psi_y[i][j][k + M*g + 3*M/4] - psi_y_p[k + M*g + 3*M/4]);
				}
			}

			alpha_x[i][j] = mult(psi_in_x_inv[i][j], b_in_x);
			alpha_y[i][j] = mult(psi_in_y_inv[i][j], b_in_y);

			for(int k = 0 ; k < MG ; k++){
				psi_x_p_out[i][j][k] = psi_x_p[k];
				psi_y_p_out[i][j][k] = psi_y_p[k];
			}

			psi_x_out = mult(psi_in_x[i][j], alpha_x[i][j]);
			psi_y_out = mult(psi_in_y[i][j], alpha_y[i][j]);
			
			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M/4 ; k++){
					// psi_x[i][j+1][k + g*M] = psi_x_out[k + g*M] + psi_x_p[k + g*M];
					// psi_x[i][j][k + g*M + M/4] = psi_x_out[k + g*M + M/4] + psi_x_p[k + g*M + M/4];
					// psi_x[i][j][k + g*M + M/2] = psi_x_out[k + g*M + M/2] + psi_x_p[k + g*M + M/2];
					psi_x[i][j+1][k + g*M + 3*M/4] = psi_x_out[k + g*M + 3*M/4] + psi_x_p[k + g*M + 3*M/4]; 

					// psi_y[i][j][k + g*M] = psi_y_out[k + g*M] + psi_y_p[k + g*M];
					// psi_y[i][j][k + g*M+M/4] = psi_y_out[k + g*M + M/4] + psi_y_p[k + g*M + M/4];
					// psi_y[i+1][j][k + g*M+M/2] = psi_y_out[k + g*M + M/2] + psi_y_p[k + g*M + M/2];
					psi_y[i+1][j][k + g*M + 3*M/4] = psi_y_out[k + g*M + 3*M/4] + psi_y_p[k + g*M + 3*M/4];
				}
			}
		}
	}

}