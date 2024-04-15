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

tuple<vector<vector <double>>, vector<vector <double>>> getLeak_psi(const int M, const int size_x, const int size_y, const int nxf, const int nyf, const int G,
								vector<vector<double>> h_x_nodes, vector<vector<double>> h_y_nodes, vector<vector<vector <double>>> psi_x, 
								vector<vector<vector <double>>> psi_y, vector<double> mi, vector<double> n, vector<double> w, vector<vector<int>> nodes_x,
								vector<vector<int>> nodes_y){

	vector<vector <double>> J_x (nyf, vector<double> (2 * G,0)); 
	vector<vector <double>> J_y (2 * G, vector<double> (nxf,0));

	vector<vector <double>> J_x_out (size_y, vector<double> (2 * G,0)); 
	vector<vector <double>> J_y_out (2 * G, vector<double> (size_x,0));

	int index = 0;

	for(int g = 0 ; g < G ; g++){
		for(int i = 0 ; i < nyf ; i++){
			for(int j = 0 ; j < M/4 ; j++){
				J_x[i][G + g] = J_x[i][G + g] + mi[j]         * psi_x[i][nxf][j + M*g]         * w[j] * h_y_nodes[i + index][0];
				J_x[i][G + g] = J_x[i][G + g] + mi[j + 3*M/4] * psi_x[i][nxf][j + 3*M/4 + M*g] * w[j + 3*M/4] * h_y_nodes[i + index][0];
			}
			J_x_out[0][G + g] = J_x_out[0][G + g] + 0.25 * J_x[i][G + g];
		}	
	}
	
	for(int g = 0 ; g < G ; g++){
		for(int i = 0 ; i < nxf ; i++){
			for(int j = 0 ; j < M/4 ; j++){
				J_y[g][i] = J_y[g][i] + n[j]         * psi_y[0][i][j + M*g]         * w[j] * h_y_nodes[0][i + index];
				J_y[g][i] = J_y[g][i] + n[j + M/4] * psi_y[0][i][j + M/4 + M*g] * w[j + M/4] * h_y_nodes[0][i + index];
			}
			J_y_out[g][0] = J_y_out[g][0] + 0.25 * J_y[g][i];
		}	
	}


	return {J_x_out, J_y_out};
}