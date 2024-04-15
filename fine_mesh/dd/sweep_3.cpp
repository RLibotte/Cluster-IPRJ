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

void sweep_3(int M, int G, int nxf, int nyf, vector<double> mi, vector<double> n, vector<double> w,
	vector<vector<vector<double>>> &psi_x_b, vector<vector<vector<double>>> &psi_y_b,
	vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y, 
	vector<vector<double>> &h_x_nodes, vector<vector<double>> &h_y_nodes,
	vector<vector<vector<double>>> &SS_x, vector<vector<vector<double>>> &SS_y, 
	vector<vector<vector<vector<double>>>> &SS_x_1, vector<vector<vector<vector<double>>>> &SS_y_1, 
	vector<vector<vector<double>>> &q_nodes, vector<vector<double>> &sigmat_z, vector<vector<vector<double>>> &sigmas1_z,
	vector<vector<int>> &zcn){

	for(int i = 0 ; i < nyf ; i++){
		for(int j = nxf - 1 ; j >= 0 ; j--){
			for(int g = 0 ; g < G ; ++g){
				for(int k = M/2 ; k < 3*M/4 ; k++){
					psi_x_b[i][j][k + g*M] = (2 * fabs(mi[k]) * psi_x[i][j+1][k + g*M]/h_x_nodes[i][j] + 2 * fabs(n[k]) * psi_y[i][j][k + g*M]/h_y_nodes[i][j] + 
										SS_x[i][j][g] + SS_x_1[i][j][g][k] + q_nodes[i][j][g]) / 
										(2 * fabs(mi[k])/h_x_nodes[i][j] + 2 * fabs(n[k])/h_y_nodes[i][j] + sigmat_z[zcn[i][j]][g]);

					psi_y_b[i][j][k + g*M] = (2 * fabs(mi[k]) * psi_x[i][j+1][k + g*M]/h_x_nodes[i][j] + 2 * fabs(n[k]) * psi_y[i][j][k + g*M]/h_y_nodes[i][j] + 
										SS_y[i][j][g] + SS_y_1[i][j][g][k] + q_nodes[i][j][g]) / 
										(2 * fabs(mi[k])/h_x_nodes[i][j] + 2 * fabs(n[k])/h_y_nodes[i][j] + sigmat_z[zcn[i][j]][g]);

					psi_x[i][j][k + g*M] = 2 * psi_x_b[i][j][k + g*M] - psi_x[i][j+1][k + g*M];
					psi_y[i+1][j][k + g*M] = 2 * psi_y_b[i][j][k + g*M] - psi_y[i][j][k + g*M];
				}
			}
		}
    }

}