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

double getdev(const int size_x, const int size_y, const int G, vector<vector<vector <double>>> &psi, vector<vector<vector <double>>> &psi_old){
	double dev = 0;
	for(int i = 0 ; i < size_y ; i++){ 
		for(int j = 0 ; j < size_x ; j++){ 
			for(int g = 0 ; g < G ; g++){
				if(fabs((psi[i][j][g]-psi_old[i][j][g])/psi[i][j][g]) > dev) dev = fabs((psi[i][j][g]-psi_old[i][j][g])/psi[i][j][g]);
			}
		}
	}
	return dev;
}