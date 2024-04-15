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

	void get_b_in(int i, int j, int G, int M, vector<double> &b_in_x, vector<double> &b_in_y,
		vector<vector<vector<double>>> &b_aux_x, vector<vector<vector<double>>> &b1_x,
		vector<vector<vector<double>>> &b_aux_y, vector<vector<vector<double>>> &b1_y,
		vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y){
	
	int index;

	for(int g = 0 ; g < G ; g++){
		for(int k = 0 ; k < M/4 ; k++){
			
			index = k + g*M;
			b_in_x[index] = psi_x[i][j][index]     - (b_aux_x[i][j][index] - b1_x[i][j][index]);
			b_in_y[index] = psi_y[i + 1][j][index] - (b_aux_y[i][j][index] - b1_y[i][j][index]);

			index += M/4;
			b_in_x[index] = psi_x[i][j + 1][index] - (b_aux_x[i][j][index] + b1_x[i][j][index]);
			b_in_y[index] = psi_y[i + 1][j][index] - (b_aux_y[i][j][index] - b1_y[i][j][index]);
	
			index += M/4;
			b_in_x[index] = psi_x[i][j + 1][index] - (b_aux_x[i][j][index] + b1_x[i][j][index]);
			b_in_y[index] = psi_y[i][j][index]     - (b_aux_y[i][j][index] + b1_y[i][j][index]);

			index += M/4;
			b_in_x[index] = psi_x[i][j][index]     - (b_aux_x[i][j][index] - b1_x[i][j][index]);
			b_in_y[index] = psi_y[i][j][index]     - (b_aux_y[i][j][index] + b1_y[i][j][index]);
		}
	}
}