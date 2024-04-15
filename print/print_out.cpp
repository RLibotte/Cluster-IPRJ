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

void print_psi_out(const int M, const int G, int wy, int wx, int nxf, int nyf, 
	vector<vector<vector<double>>> psi_x, vector<vector<vector<double>>> psi_y){

	double temp = 0;
	for(int g = 0 ; g < G ; g++){
		temp = 0;
		for(int i = 0 ; i < wy ; i++){
			for(int j = 0 ; j < M/4 ; j++){
				temp = temp + psi_x[i][nxf][j + g*M] + psi_x[i][nxf][j + 3*M/4 + g*M];
			}
		}
		printf("face leste (g = %d) = %.10f \n\n", g, temp/(wy * (M/4)));
	}
	
	for(int g = 0 ; g < G ; g++){
		temp = 0;
		for(int i = 0 ; i < wx ; i++){
			for(int j = 0 ; j < M/4 ; j++){
				temp = temp + psi_y[0][i][j + g*M] + psi_y[0][i][j + M/4 + g*M];
			}
		}
		printf("face norte = %.10f \n\n", temp/(wx * (M/4)));
	}
	
 	for(int g = 0 ; g < G ; g++){
 		temp = 0;
 		for(int i = 0 ; i < wy ; i++){
			for(int j = M/4 ; j < M/2 ; j++){
				temp = temp + psi_x[i][0][j + g*M] + psi_x[i][0][j + M/4 + g*M];
			}
		}
		printf("face oeste = %.10f \n\n", temp/(wy * (M/4)));
 	}
	
	for(int g = 0 ; g < G ; g++){
		temp = 0;
		for(int i = 0 ; i < wx ; i++){
			for(int j = M/2 ; j < 3*M/4 ; j++){
				temp = temp + psi_y[nyf][i][j + g*M] + psi_y[nyf][i][j + M/4 + g*M];
			}
		}
		printf("face sul   = %.10f \n\n", temp/(wx * (M/4)));

	}
}

void print_leakage_out(const int M, const int G, int x_size, int y_size, 
	vector<vector<double>> J_x, vector<vector<double>> J_y){

	cout << "Fuga X" << endl << endl;
	for(int g = 0 ; g < G ; g++){
		cout << "Grupo " << g+1 << endl << endl;
		for(int j = 0 ; j < y_size ; j++){
			printf("%.8e \t", J_x[j][g]);
			printf("%.8e \n", J_x[j][G + g]);
		}
		cout << endl;
	}
	
	cout << endl;

	cout << "Fuga Y" << endl << endl;;
	for(int g = 0 ; g < G ; g++){
		cout << "Grupo " << g+1 << endl << endl;
		for(int j = 0 ; j < x_size ; j++){
			printf("%.8e \t", J_y[g][j]);
		}
		cout << endl;
		for(int j = 0 ; j < x_size ; j++){
			printf("%.8e \t", J_y[g + G][j]);
		}
		cout << endl << endl;
	}
	cout << endl;

	double temp_f = 0;

	for(int g = 0 ; g < G ; g++){
		temp_f = 0;
		for(int i = 0 ; i < y_size ; i++){
			temp_f = temp_f +  J_x[i][g + G];
		}
		cout << "grupo " << g+1 << endl;
		printf("Fuga direita = %.8e \n", temp_f);

	}
	
	cout << endl;

	for(int g = 0 ; g < G ; g++){
		temp_f = 0;
		for(int i = 0 ; i < x_size ; i++){
			temp_f = temp_f +  J_y[g][i];
		}
		cout << "grupo " << g+1 << endl;
		printf("Fuga superior = %.8e \n", temp_f);
	}
	cout << endl;
}