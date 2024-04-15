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

tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<vector<double>>>,	
	  vector<vector<vector<double>>>, vector<vector<vector<double>>>> 

	get_eig(int M, int G, int Z, int x_size, int y_size, int nxf, int nyf, vector<vector<int>> &zone_config,
	vector<double> &mi, vector<double> &n, vector<double> &w, vector<vector<double>> &sigmat_z, 
	vector<vector<vector<double>>> &sigmas0_z, vector<vector<vector<double>>> &sigmas1_z, 
	vector<vector<int>> &nodes_x, vector<vector<int>> &nodes_y){

	const int MG = M*G;

	double index_x = 0, index_y = 0;

	vector<vector<vector<double>>> A_x (Z, vector<vector<double>> (MG, vector<double> (MG)));
	vector<vector<vector<double>>> A_y (Z, vector<vector<double>> (MG, vector<double> (MG)));
	vector<vector<vector<double>>> A_p_inv_z (Z, vector<vector<double>> (MG, vector<double> (MG)));

	vector<vector<vector<vector<double>>>> A_p_inv (y_size, vector<vector<vector<double>>> (x_size, vector<vector<double>> (MG, vector<double>(MG))));
	vector<vector<vector<vector<double>>>> A_p_inv_out (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>> (MG, vector<double>(MG))));

	vector<vector<double>> D_x (Z, vector <double>(MG));
	vector<vector<double>> D_y (Z, vector <double>(MG));
	
	vector<vector<vector<double>>> V_x (Z, vector<vector<double>> (MG, vector <double>(MG)));
	vector<vector<vector<double>>> V_y (Z, vector<vector<double>> (MG, vector <double>(MG)));

	int aux;


	for(int l = 0 ; l < Z ; l++){
		for(int g = 0 ; g < G ; ++g){	
			for(int i = 0 ; i < M ; ++i){
				aux = i + g*M;
				for(int h = 0 ; h < G ; ++h){
					for(int j = 0 ; j < M ; ++j){
						A_p_inv_z[l][aux][j + h*M] = - 0.25 * w[j] * (sigmas0_z[l][h][g] + 3.0 * sigmas1_z[l][h][g] * (mi[i] * mi[j] + n[i] * n[j]));
					}
				}

				A_p_inv_z[l][aux][aux] = sigmat_z[l][g] + A_p_inv_z[l][aux][aux];

				for(int j = 0 ; j < M*G ; j++){
					A_x[l][aux][j] = A_p_inv_z[l][aux][j] / mi[i];
					A_y[l][aux][j] = A_p_inv_z[l][aux][j] / n[i];
				}	
			}
		}

		auto[V_x_temp, D_x_temp] = eigen(A_x[l], MG);
		auto[V_y_temp, D_y_temp] = eigen(A_y[l], MG);

		V_x[l] = V_x_temp;
		V_y[l] = V_y_temp;

		for(int i = 0 ; i < MG ; i++){
			D_x[l][i] = 1.0/D_x_temp[i];
			D_y[l][i] = 1.0/D_y_temp[i];
		}

		A_p_inv_z[l] = inverse(A_p_inv_z[l], MG);
	}

	return {D_x, D_y, V_x, V_y, A_p_inv_z};
	
}