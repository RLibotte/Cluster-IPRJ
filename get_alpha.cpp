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

vector<double> get_alpha(int N, int G, int k, 
    vector<vector<double>> psi, vector<vector<double>> psi_p, vector<vector<double>> arr){

    vector<double> psi_temp (N*G);
    
    for (int i = 0 ; i < G ; i++){
        for(int j = 0 ; j < N/2 ; j++){
             psi_temp[j + i*N] = psi[j + i*N][k] - psi_p[k][j + i*N];
        }
        for(int j = N/2 ; j < N ; j++){
            psi_temp[j + i*N] = psi[j + i*N][k + 1] - psi_p[k][j + i*N];
        }
    }

    return mult(arr, psi_temp);
}