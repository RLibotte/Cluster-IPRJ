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

vector<vector<double>> get_SS(int N, int G, int k, vector<vector<vector<double>>> sigmas0, 
    vector<vector<vector<double>>> sigmas1,	vector<double> mi, vector<double> w, vector<double> psi_m, vector<vector<double>> SS){

    double sum, sum1;

    for(int i = 0 ; i < G ; i++){
        sum = 0;
        sum1 = 0;
        #pragma omp parallel for collapse (2) reduction(+:sum)
        for(int j = 0 ; j < G ; j++){
            for(int l = 0 ; l < N ; l++){
                sum += (sigmas0[k][j][i]/2)*(w[l] * psi_m[l + j*N]);
                sum1 += (sigmas1[k][j][i])*(w[l] * psi_m[l + j*N] * mi[l]);
            }
        }
        SS[i][0] = sum;
        SS[i][1] = sum1;
    }

    return SS;

}