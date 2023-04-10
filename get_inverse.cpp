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

void get_inverse(const int NG, int R, std::vector<std::vector<std::vector<double>>> arr, std::vector<std::vector<std::vector<double>>> &Am_inv){

	Am_inv.resize(R, std::vector<std::vector<double>> (NG, std::vector<double> (NG)));

	for(int i = 0 ; i < R ; i++){
		Am_inv[i] = inverse(arr[i], NG);
	}
}