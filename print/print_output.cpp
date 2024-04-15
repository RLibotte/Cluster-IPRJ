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

void print_output(int N, int x_size, int y_size, int nxf, int nyf, double hx, double hy, int G, vector<vector<vector<double>>> scalar_flux, vector<vector<vector<double>>> t_abs,
				  vector<vector<vector<double>>> scalar_flux_out, vector<vector<vector<double>>> t_abs_out, string method, string model){

	string filename = "../Output/Model_" + model + "_S" + to_string(N) + "_" + method + "_scalar_flux_out_" + to_string(nxf) + "x" + to_string(nyf) + ".txt";

	std::ofstream scalar_flux_out_f (filename);
	for(int g = 0 ; g < G ; g++){
		scalar_flux_out_f << "G = " << g+1 << endl << endl;
		for(int i = y_size-1 ; i >= 0 ; i--){
			for(int j = 0 ; j < x_size ; j++){
				scalar_flux_out_f << std::scientific << scalar_flux_out[i][j][g] << " \t";
			}
			scalar_flux_out_f << endl;
		}
		scalar_flux_out_f << endl;
	}
	scalar_flux_out_f.close();

	//////////////////////////////////////////////////////////////////////////////////////////////

	filename = "../Output/Model_" + model + "_S" + to_string(N) + "_" + method + "_t_abs_out_" + to_string(nxf) + "x" + to_string(nyf) + ".txt";;

	std::ofstream t_abs_out_f (filename);
	for(int g = 0 ; g < G ; g++){
		t_abs_out_f << "G = " << g+1 << endl << endl;
		for(int i = 0 ; i < y_size ; i++){
			for(int j = 0 ; j < x_size ; j++){
				t_abs_out_f << std::scientific << t_abs_out[i][j][g] << " \t";
			}
			t_abs_out_f << endl;
		}
		t_abs_out_f << endl;
	}
	t_abs_out_f.close();

	//////////////////////////////////////////////////////////////////////////////////////////////

	filename = "../Output/Model_" + model + "_S" + to_string(N) + "_" + method + "_scalar_flux_" + to_string(nxf) + "x" + to_string(nyf) + ".txt";;

	std::ofstream scalar_flux_f (filename);
	for(int g = 0 ; g < G ; g++){
		scalar_flux_f << "G = " << g+1 << endl << endl;
		for(int i = nyf-1 ; i >= 0 ; i--){
			for(int j = 0 ; j < nxf ; j++){
				scalar_flux_f << std::scientific << scalar_flux[i][j][g] << " \t";
			}
			scalar_flux_f << endl;
		}
		scalar_flux_f << endl;
	}
	scalar_flux_f.close();

	//////////////////////////////////////////////////////////////////////////////////////////////

	filename = "../Output/Model_" + model + "_S" + to_string(N) + "_" + method + "_t_abs_" + to_string(nxf) + "x" + to_string(nyf) + ".txt";;

	std::ofstream t_abs_f (filename);
	for(int g = 0 ; g < G ; g++){
		t_abs_f << "G = " << g+1 << endl << endl;
		for(int i = 0 ; i < nyf ; i++){
			for(int j = 0 ; j < nxf ; j++){
				t_abs_f << std::scientific << t_abs[i][j][g] << " \t";
			}
			t_abs_f << endl;
		}
		t_abs_f << endl;
	}
	t_abs_f.close();

	//////////////////////////////////////////////////////////////////////////////////////////////

	filename = "../Output/" + method + "_scalar_flux_chart.txt";

	std::ofstream scalar_flux_chart (filename);
	scalar_flux_chart << nxf << " \t" << nyf <<  " \t" << hx <<  " \t" << hy <<  " \t" << G << endl << endl;
	for(int g = 0 ; g < G ; g++){
		for(int i = nyf-1 ; i >= 0 ; i--){
			for(int j = 0 ; j < nxf ; j++){
				scalar_flux_chart << std::scientific << scalar_flux[i][j][g] << " \t";
			}
			scalar_flux_chart << endl;
		}
		scalar_flux_chart << endl;
	}
	scalar_flux_chart.close();

	// printf("Output files written succesfully.\n");
}