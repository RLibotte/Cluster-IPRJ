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

void get_psi_l(int N, int G, int NG, int k, vector<vector<double>> &psi, vector<vector<double>> ni,
	vector<vector<double>> temp, vector<vector<vector<double>>> V, vector<double> h, 
	vector<vector<double>> psi_p, vector<double> alpha, vector<double> psi_temp){

	for(int i = 0 ; i < NG ; i++){
        for(int j = 0 ; j < NG ; j++){
            if(ni[j][k] > 0){
                temp[i][j] = V[k][i][j];
            } else {
                temp[i][j] = V[k][i][j] * exp(h[k] / ni[j][k]);
            }
        }
    }

    psi_temp = mult(temp, alpha);

    for(int i = 0 ; i < G ; i++){
        for(int j = N/2 ; j < N ; j++){
            psi[j + i*N][k] = psi_temp[j + i*N] + psi_p[k][j + i*N];
        }
    }
}

void get_psi_r(int N, int G, int NG, int k, vector<vector<double>> &psi, vector<vector<double>> ni,
	vector<vector<double>> temp, vector<vector<vector<double>>> V, vector<double> h, 
	vector<vector<double>> psi_p, vector<double> alpha, vector<double> psi_temp){

	for(int i = 0 ; i < NG ; i++){
	    for(int j = 0 ; j < NG ; j++){
	        if(ni[j][k] > 0){
	            temp[i][j] = V[k][i][j] * exp(-h[k] / ni[j][k]);
	        } else {
	            temp[i][j] = V[k][i][j];
	        }
	    }
	}

	psi_temp = mult(temp, alpha);

	for(int i = 0 ; i < G ; i++){
	    for(int j = 0 ; j < N/2 ; j++){
	        psi[j + i*N][k+1] = psi_temp[j + i*N] + psi_p[k][j + i*N];
	    }
	}
}