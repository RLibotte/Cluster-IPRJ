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

void print_output(int G, int nxf, vector<vector<double>> scalar_flux_out, vector<double> h){

	string filename = "../output/scalar_flux.txt";
	double pos = 0;

	ofstream scalar_flux (filename);
	scalar_flux << scientific << G << "," << endl;
	scalar_flux << scientific << nxf+1 << "," << endl;

	for(int i = 0 ; i < nxf ; i++){
		scalar_flux << std::scientific << h[i] << ", \t";
	}
	scalar_flux << std::scientific << h[nxf] << endl;

	for(int g = 0 ; g < G ; g++){
		for(int i = 0 ; i < nxf ; i++){
				scalar_flux << std::scientific << scalar_flux_out[g][i] << ", \t";
			}
			scalar_flux << std::scientific << scalar_flux_out[g][nxf] << "\t";
			scalar_flux << endl;
		}
	scalar_flux.close();
}

vector<double> get_pos_arr(int R, int nxf, vector<double> h, vector<int> nodes){

	double total = 0;
    double pos = 0;

    vector<double> pos_arr (nxf+1);
    vector<double> pos_r (nxf+1);

    int ind = 1;
    for(int i = 0 ; i < R ; i++){
        for(int j = 0 ; j < nodes[i] ; j++){
            pos_arr[ind] = pos_arr[ind - 1] + h[i]/nodes[i];
            ind++;
        }
    }

    double pos_temp = 0;
    ind = 1;
    for(int i = 0 ; i < R; i++){
        for(int j = 0 ; j < nodes[i] - 1; j++){
            pos_r[ind] = pos_temp + h[i]/nodes[i];
            pos_temp += h[i]/nodes[i];
            ind++;
        }
        pos_temp = 0;
        ind++;
    }
    pos_r[ind - 1] = pos_r[ind - 2] + h[R-1]/nodes[R-1];

    return pos_arr;
}