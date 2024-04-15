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

void get_transverse_leakage(int nxf, int nyf, int M, int G, int lb, int ub,
	vector<vector<vector<double>>>  &psi_x, vector<vector<vector<double>>>  &psi_y,
	vector<double> &mi, vector<double> &n, vector<vector<double>> &h_x, vector<vector<double>> &h_y,
	vector<vector<vector<double>>> &l_x, vector<vector<vector<double>>> &l_y){

	int aux, i, j ,g ,k;

	for(i = 0 ; i < nyf ; i++){
		for(j = 0 ; j < nxf ; j++){
			for(g = 0 ; g < G ; g++){
				for(k = lb ; k < ub ; k++){
					aux = k + g*M;

					l_x[i][j][aux]  = (mi[k]/h_x[i][j]) * (psi_x[i][j+1][aux] - psi_x[i][j][aux]);
					l_y[i][j][aux]  =  (n[k]/h_y[i][j]) * (psi_y[i][j][aux] - psi_y[i+1][j][aux]);
				}
			}
		}
	}
}