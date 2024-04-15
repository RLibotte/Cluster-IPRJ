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

void get_psi_p(int nxf, int nyf, int G, int M, vector<vector<vector<double>>> &l_x, vector<vector<vector<double>>> &l_y,
	vector<double> &l_x_temp, vector<double> &l_y_temp, vector<vector<vector<double>>> &q, vector<double> &psi_x_p,
	vector<double> &psi_y_p, vector<vector<vector<vector<double>>>> &A_p_lu, vector<double> &b_in_x,vector<double> &b_in_y,
	vector<vector<vector<double>>> &alpha_x, vector<vector<vector<double>>> &alpha_y, vector<vector<vector<vector<double>>>> &psi_in_x_inv,
	vector<vector<vector<vector<double>>>> &psi_in_y_inv, vector<vector<vector<double>>> &psi_x_p_out, vector<vector<vector<double>>> &psi_y_p_out,
	vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y){

	int index;

	for(int i = 0 ; i < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; ++g){
				for(int k = 0 ; k < M ; k++){
					index = k + g*M;

					l_x_temp[index] = q[i][j][g] - l_x[i][j][index];
					l_y_temp[index] = q[i][j][g] - l_y[i][j][index];
				}
			}

			psi_x_p = mult(A_p_lu[i][j], l_y_temp);			
			psi_y_p = mult(A_p_lu[i][j], l_x_temp);

			for(int g = 0 ; g < G ; g++){
				for(int k = 0 ; k < M/4 ; k++){
					index = k + g*M;
					b_in_x[index] = psi_x[i][j][index] - psi_x_p[index];
					b_in_y[index] = psi_y[i+1][j][index] - psi_y_p[index];

					index += M/4;
					b_in_x[index] = psi_x[i][j + 1][index] - psi_x_p[index];
					b_in_y[index] = psi_y[i+1][j][index] - psi_y_p[index];

					index += M/4;
					b_in_x[index] = psi_x[i][j + 1][index] - psi_x_p[index];
					b_in_y[index] = psi_y[i][j][index] - psi_y_p[index];

					index += M/4;
					b_in_x[index] = psi_x[i][j][index] - psi_x_p[index];
					b_in_y[index] = psi_y[i][j][index] - psi_y_p[index];
				}
			}

			alpha_x[i][j] = mult(psi_in_x_inv[i][j], b_in_x);	
			alpha_y[i][j] = mult(psi_in_y_inv[i][j], b_in_y);	

			psi_x_p_out[i][j] = psi_x_p;
			psi_y_p_out[i][j] = psi_y_p;
		}
	}
}