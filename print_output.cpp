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

void print_output(int G, int nxf, vector<vector<double>> scalar_flux_out, vector<double> h){

	string filename = "../output/scalar_flux_plot.txt";
	double pos = 0;

	std::ofstream scalar_flux (filename);
	scalar_flux << std::scientific << G << "," << endl;
	scalar_flux << std::scientific << nxf+1 << "," << endl;

	for(int i = 0 ; i < nxf ; i++){
		scalar_flux << std::scientific << h[i] << ", \t";
	}
	scalar_flux << std::scientific << h[nxf] << endl;

	// for(int i = 0 ; i < nxf ; i++){
	// 	scalar_flux << std::scientific << pos << ", \t";
	// 	pos += h[i];
	// }
	// scalar_flux << std::scientific << pos << endl;


	for(int g = 0 ; g < G ; g++){
		for(int i = 0 ; i < nxf ; i++){
				scalar_flux << std::scientific << scalar_flux_out[g][i] << ", \t";
			}
			scalar_flux << std::scientific << scalar_flux_out[g][nxf] << "\t";
			scalar_flux << endl;
		}
	scalar_flux.close();
}