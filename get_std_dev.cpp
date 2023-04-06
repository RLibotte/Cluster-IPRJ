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

void get_std_dev(int iter, double sum, vector<double> avg){

	double std_dev = 0;

    for(int i = 0 ; i < iter ; i++){
        std_dev += pow(avg[i] - sum/iter, 2);
    }

    std_dev = sqrt(std_dev/iter);
    printf("Standard Deviation = %.6f s\n", std_dev);
}