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

vector<vector <double>> get_phi(vector<vector<double>> phi, vector<vector<double>> psi, vector<double> w, int nxf, int G, int N){
	for(int i = 0 ; i < G ; i++){
        for(int j = 0 ; j < nxf + 1 ; j++){
            phi[i][j] = 0;
        }
    }

    for(int i = 0 ; i < G ; i++){
        for(int j = 0 ; j < nxf + 1 ; j++){
            for(int k = 0 ; k < N ; k++){
                phi[i][j] = phi[i][j] + psi[i*N + k][j]*w[k];
            }
            phi[i][j] = phi[i][j]/2;   
        }
    }

    return phi;
}