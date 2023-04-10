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

vector<double> get_psi_m(int N, int G, int k, vector<vector<double>> ni, 
	vector<vector<vector<double>>> V, vector<double> alpha, vector<double> h, vector<double> psi_p){

 	vector<double> psi_m (N*G, 0);
 
    double sum;

	for(int i = 0 ; i < N * G ; i++){
        sum = 0;
        #pragma omp parallel for reduction (+:sum)
        for(int j = 0 ; j < N * G; j++){
            if(ni[j][k] > 0){
                sum += (alpha[j] * V[k][i][j] * ni[j][k] * (exp(-h[k]/ni[j][k]) - 1));
            } else {
                sum += (alpha[j] * V[k][i][j] * ni[j][k] * (1 - exp(h[k]/ni[j][k])));
            }
        }
        psi_m[i] = -(1 / h[k]) * sum + psi_p[i];
    }
    return psi_m;
}

vector<double> get_psi_m_h(int N, int G, int k, vector<vector<double>> ni, 
    vector<vector<vector<double>>> V, vector<double> alpha, vector<double> h, vector<double> psi_p){

    vector<double> psi_m (N*G, 0);
 
    double sum;

    for(int i = 0 ; i < (N*G)/2 ; i++){
        sum = 0;
        #pragma omp parallel for reduction (+:sum)
        for(int j = 0 ; j < (N*G)/2; j++){
            sum += (alpha[j] * V[k][i][j] * ni[j][k] * (1 - exp(-h[k]/ni[j][k])));
            sum += (alpha[j + N/2] * V[k][i + N/2][j] * ni[j][k] * (1 - exp(-h[k]/ni[j][k])));
        }
        psi_m[i] = -(1 / h[k]) * sum + psi_p[i];
    }

    for(int i = 0 ; i < (N*G)/2 ; i++){
        sum = 0;
        #pragma omp parallel for reduction (+:sum)
        for(int j = 0 ; j < (N*G)/2; j++){
            sum += (alpha[j] * V[k][i+ N/2][j] * ni[j][k] * (exp(-h[k]/ni[j][k])) - 1);
            sum += (alpha[j + N/2] * V[k][i][j] * ni[j][k] * (exp(-h[k]/ni[j][k])) - 1);
        }
        psi_m[i+N/2] = -(1 / h[k]) * sum + psi_p[i];
    }
    return psi_m;
}