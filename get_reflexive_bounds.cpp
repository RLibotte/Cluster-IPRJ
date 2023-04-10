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

std::vector<std::vector<double>> reflexive_l(std::vector<std::vector<double>> psi, int N, int G){

    for(int i = 0 ; i < G ; i++){
        for(int j = 0 ; j < N/2 ; j++){
            psi[i*N + j][0] = psi[i*N + j + N/2][0];
        }
    }

    return psi;

}

std::vector<std::vector<double>> reflexive_r(std::vector<std::vector<double>> psi, int N, int G, int nxf){

    for(int i = 0 ; i < G ; i++){
        for(int j = 0 ; j < N/2 ; j++){
            psi[i*N + j + N/2][nxf] = psi[i*N + j][nxf];
        }
    }

    return psi;

}