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

using namespace std;

double get_dev(vector<vector<double>> phi, vector<vector<double>> phi_old, int G, int nxf){

	double new_dev = 0;

	#pragma omp parallel for shared(new_dev)
	for(int i = 0 ; i < G ; i++){
		for(int j = 0 ; j < nxf+1 ; j++){
			if(fabs((phi[i][j] - phi_old[i][j])/phi_old[i][j]) > new_dev) new_dev =  fabs((phi[i][j] - phi_old[i][j])/phi_old[i][j]);
		}
	}

	for(int i = 0 ; i < G ; i++){
        for(int j = 0 ; j < nxf+1 ; j++){
            phi_old[i][j] = phi[i][j];
        }
    }

	return new_dev;

}