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

void get_scalar_flux(const int &M, const int &G, int &nxf, int &nyf, vector<vector<vector<double>>> &scalar_flux_x,
					 vector<vector<vector<double>>> &scalar_flux_x_out, vector<vector<vector<double>>> &psi_b_x,
					 vector<vector<int>> &nodes_x, vector<vector<int>> &nodes_y, vector<double> &w, int x_size, int y_size){

	int index_x = 0, index_y = 0;

	for(int i = 0 ; i < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; ++g){
				scalar_flux_x[i][j][g] = 0;
				for(int k = 0 ; k < M ; k++){
					scalar_flux_x[i][j][g] = scalar_flux_x[i][j][g] + psi_b_x[i][j][k + g*M] * w[k];
				}
				scalar_flux_x[i][j][g] = 0.25 * scalar_flux_x[i][j][g];
			}
		}
	}

	for(int i = 0 ; i < y_size ; i++){
		for(int j = 0 ; j < x_size ; j++){
			for(int g = 0 ; g < G ; ++g){
				scalar_flux_x_out[i][j][g] = 0;
				for(int k = 0 ; k < nodes_y[i][j] ; k++){
					for(int l = 0 ; l < nodes_x[i][j] ; l++){
						scalar_flux_x_out[i][j][g] = scalar_flux_x_out[i][j][g] + scalar_flux_x[k + index_y][l + index_x][g]; 
					}
				}
				scalar_flux_x_out[i][j][g] = scalar_flux_x_out[i][j][g] / (nodes_x[i][j] * nodes_y[i][j]);
			}
			index_x = index_x + nodes_x[i][j];
		}
		index_x = 0;
		index_y = index_y + nodes_y[i][0];
	}	
	index_y = 0;
}