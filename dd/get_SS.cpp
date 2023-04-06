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

std::vector<std::vector<std::vector<double>>> get_SS(std::vector<std::vector<std::vector<double>>> SS, std::vector<std::vector<std::vector <double>>> sigmas0, 
            std::vector<std::vector<std::vector <double>>> sigmas1, std::vector<double> mi, std::vector<double> w, 
            std::vector<std::vector<double>> psi_m, int nxf, int G, int N, int L){

	// #pragma omp parallel for collapse(2) shared(SS)
    for(int k = 0 ; k < nxf ; k++){
        for(int i = 0 ; i < G ; i++){
            for(int j = 0 ; j < L ; j++){
                SS[i][k][j] = 0;
            }
        }
    }
        
    // #pragma omp parallel for shared (SS) collapse(4)
    for(int k = 0 ; k < nxf ; k++){
        for(int i = 0 ; i < G ; i++){
            for(int j = 0 ; j < G ; j++){
                for(int m = 0 ; m < N ; m++){
                    SS[i][k][0] += (sigmas0[k][j][i]/2) * (w[m] * psi_m[m + j * N][k]);
                    SS[i][k][1] += (sigmas1[k][j][i]) * (w[m] * psi_m[m + j * N][k] * mi[m]);
                }
            }
        }
    }

    return SS;
}