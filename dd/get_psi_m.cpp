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

std::vector<std::vector<double>> get_psi_m(std::vector<std::vector<double>> psi, std::vector<std::vector<double>> psi_m, int NG, int nxf){

   	#pragma omp parallel for collapse(2)
    for(int i = 0 ; i < NG ; i++){
        for(int j = 0 ; j < nxf ; j++){
           psi_m[i][j] = (psi[i][j] + psi[i][j+1])/2;
        }
    }

    return psi_m;

}