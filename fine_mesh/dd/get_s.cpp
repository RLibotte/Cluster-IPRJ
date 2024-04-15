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

void get_s(const int M, const int G, const int nxf, const int nyf, 
			vector<double> &w, vector<double> &mi, vector<double> &n, 
			vector<vector<vector<double>>> &sigmas0_z, vector<vector<vector<double>>> &sigmas1_z,
			vector<vector<vector<double>>> &psi_x_b, vector<vector<vector <double>>> &psi_y_b,
			vector<vector<vector<double>>> &SS_x, vector<vector<vector<double>>> &SS_y,
			vector<vector<vector<vector<double>>>> &SS_x_1, vector<vector<vector<vector<double>>>> &SS_y_1,
			vector<vector<int>> &zcn){

	double temp_x = 0, temp_y = 0;

	for(int i = nyf - 1 ; i >= 0 ; i--){
		for(int j = 0 ; j < nxf ; j++){
			for(int g = 0 ; g < G ; ++g){
	    		SS_x[i][j][g] = 0;
				SS_y[i][j][g] = 0;

	            for(int k = 0 ; k < M ; k++){
					SS_x_1[i][j][g][k] = 0;
					SS_y_1[i][j][g][k] = 0;
				}

	    		for(int h = 0 ; h < G ; h++){
	        		for(int l = 0 ; l < M ; l++){
		            	SS_x[i][j][g] = SS_x[i][j][g] + 0.25 * sigmas0_z[zcn[i][j]][h][g] * w[l] * psi_x_b[i][j][l + h*M];
		            	SS_y[i][j][g] = SS_y[i][j][g] + 0.25 * sigmas0_z[zcn[i][j]][h][g] * w[l] * psi_y_b[i][j][l + h*M];
		            }

		    		for(int k = 0 ; k < M/2 ; k++){		
					temp_x = 0;
					temp_y = 0;

	        			for(int l = 0 ; l < M ; l++){
	        				temp_x += w[l] * psi_x_b[i][j][l + h*M] * (mi[k] * mi[l] + n[k] * n[l]);
			            	temp_y += w[l] * psi_y_b[i][j][l + h*M] * (mi[k] * mi[l] + n[k] * n[l]);
		    			}

		    			SS_x_1[i][j][g][k] += temp_x * 0.75 * sigmas1_z[zcn[i][j]][h][g];
		    			SS_y_1[i][j][g][k] += temp_y * 0.75 * sigmas1_z[zcn[i][j]][h][g];

		    			SS_x_1[i][j][g][k + M/2] -= temp_x * 0.75 * sigmas1_z[zcn[i][j]][h][g]; /////////////////////////////
		    			SS_y_1[i][j][g][k + M/2] -= temp_y * 0.75 * sigmas1_z[zcn[i][j]][h][g]; /////////////////////////////
		            }
		        }
		    }
		}
	}
}