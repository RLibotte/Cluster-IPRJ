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

tuple<int, int, 
	vector<vector<vector<double>>>, 
	vector<vector<double>>,	
	vector<vector<double>>, 
	vector<vector<int>>,
	int, int, double, const int> 
	
	get_arrays(const int M, const int G, const int x_size, const int y_size, 
	vector<vector<vector<double>>> &q_n, vector<vector<int>> &nodes_x, vector<vector<int>> &nodes_y,
	vector<vector<double>> &h_x_n, vector<vector<double>> &h_y_n, vector<vector<int>> &zone_config){

	int nxf = 0, nyf = 0, z = 0, index;
	double dev_x = 1;

	const int MG = M*G;

	for(int i = 0 ; i < x_size ; ++i){
		nxf = nxf + nodes_x[0][i];
	}

	for(int i = 0 ; i < y_size ; ++i){
		nyf = nyf + nodes_y[i][0];
	}

	vector<vector<vector<double>>> q (nyf, vector<vector<double>> (nxf, vector <double>(G)));

	vector<vector <double>> h_x(nyf , vector<double> (nxf)); 
	vector<vector <double>> h_y(nyf , vector<double> (nxf)); 

	vector<vector <int>> zcn(nyf , vector<int> (nxf)); 

	int index_x = 0, index_y = 0;

	for(int i = 0 ; i < x_size ; ++i){
		for(int j = 0 ; j < y_size ; ++j){
			for(int k = 0 ; k < nodes_y[j][i] ; ++k){
				for(int l = 0 ; l < nodes_x[j][i] ; ++l){
					for(int g = 0 ; g < G ; ++g){
						q[k + index_y][l + index_x][g] = q_n[j][i][g];
						h_x[k + index_y][l + index_x] = h_x_n[j][i]/nodes_x[j][i];
						h_y[k + index_y][l + index_x] = h_y_n[j][i]/nodes_y[j][i];
					}
					zcn[k + index_y][l + index_x] = zone_config[j][i];
				}
			}
			index_y += nodes_y[j][i];
		}
		index_y = 0;
		index_x += nodes_x[0][i];
	}


	return {nxf, nyf, q, h_x, h_y, zcn, z, index, dev_x, MG};
}