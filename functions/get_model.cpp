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

tuple<const int, const int, const int, int, int, int, int, vector<vector <int>>, vector<vector <int>>, 
	  vector<vector<double>>, vector<vector<vector <double>>>,
      vector<vector<vector <double>>>, vector<vector<vector <double>>>, 
      vector<vector <double>>, vector<vector <double>>, const int, vector<vector<int>>>  
      get_model(string model, const int wx, const int wy){

	int x_size, y_size, G, Z, ref_o, ref_s, ref_l, ref_n;

	/////////////////////////////////////////////////////////////////////////////////////

	vector<double> data;
	string filename = "../../models_2d/model_" + model + ".txt";
	std::ifstream ifile(filename, std::ios::in);
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input Model file!\n";
    }
    double num = 0.0;
    while (ifile >> num) {
    	data.push_back(num);
    }
    
    /////////////////////////////////////////////////////////////////////////////////////

    x_size = int(data[0]);
    y_size = int(data[1]);
    
    G = int(data[2]);
    Z = int(data[3]);

    ref_o = int(data[4]);
    ref_s = int(data[5]);
    ref_l = int(data[6]);
    ref_n = int(data[7]);

    vector<vector<int>> zone_config (y_size, vector<int> (x_size));

	vector<vector<int>> nodes_x (y_size, vector<int> (x_size));
	vector<vector<int>> nodes_y (y_size, vector<int> (x_size));

	vector<vector<double>> sigmat_z (Z, vector<double> (G));

	vector<vector<vector<double>>> sigmas0_z (Z, vector<vector<double>> (G, vector<double> (G)));
	vector<vector<vector<double>>> sigmas1_z (Z, vector<vector<double>> (G, vector<double> (G)));

	vector<vector<double>> h_x (y_size, vector<double> (x_size));
	vector<vector<double>> h_y (y_size, vector<double> (x_size));

	vector<vector<vector<double>>> q (y_size, vector<vector<double>> (x_size, vector<double> (G)));

	int a = 8;

    for(int i = 0 ; i < y_size ; i++){
    	for(int j = 0 ; j < x_size ; j++){
    		zone_config[i][j] = int(data[a]) - 1;
    		a++;
    	}
    }

    for(int i = 0 ; i < y_size ; i++){
    	for(int j = 0 ; j < x_size ; j++){
    		nodes_x[i][j] = int(data[a]) * wx;
    		a++;
    	}
    }

    for(int i = 0 ; i < y_size ; i++){
    	for(int j = 0 ; j < x_size ; j++){
    		nodes_y[i][j] = int(data[a]) * wy;
    		a++;
    	}
    }

    for(int i = 0 ; i < Z ; i++){
    	for(int j = 0 ; j < G ; j++){
    		sigmat_z[i][j] = data[a];
    		a++;
    	}
    }

    for(int i = 0 ; i < Z ; i++){
    	for(int j = 0 ; j < G ; j++){
    		for(int k = 0 ; k < G ; k++){
    			sigmas0_z[i][k][j] = data[a]; /////////////////////////////////////////////
    			a++;
    		}
    	}
    }

    for(int i = 0 ; i < Z ; i++){
    	for(int j = 0 ; j < G ; j++){
    		for(int k = 0 ; k < G ; k++){
    			sigmas1_z[i][k][j] = data[a]; ///////////////////////////////////////////////
    			a++;
    		}
    	}
    }

    for(int i = 0 ; i < y_size ; i++){
    	for(int j = 0 ; j < x_size ; j++){
    		h_x[i][j] = data[a];
    		a++;
    	}
    }

    for(int i = 0 ; i < y_size ; i++){
    	for(int j = 0 ; j < x_size ; j++){
    		h_y[i][j] = data[a];
    		a++;
    	}
    }

    for(int k = 0 ; k < G ; k++){
    	for(int i = 0 ; i < y_size ; i++){
	    	for(int j = 0 ; j < x_size ; j++){
	    		q[i][j][k] = data[a];
	    		a++;
	    	}
    	}
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////

	// vector<vector<vector <double>>> sigmat_n (y_size, vector<vector<double>> (x_size, vector<double> (G)));
	// vector<vector<vector<vector<double>>>> sigmas0_n (y_size, vector<vector<vector<double>>> (x_size, vector<vector<double>>(G, vector<double> (G))));
	// vector<vector<vector<vector<double>>>> sigmas1_n (y_size, vector<vector<vector<double>>> (x_size, vector<vector<double>>(G, vector<double> (G))));

	// for(int i = 0 ; i < y_size; i++){
	// 	for(int j = 0 ; j < x_size ; j++){
	// 		for(int k = 0 ; k < G ; k++){
	// 			sigmat_n[i][j][k] = sigmat_z[zone_config[i][j]][k];
	// 			for(int l = 0 ; l < G ; l++){
	// 				sigmas0_n[i][j][k][l] = sigmas0_z[zone_config[i][j]][k][l];
	// 				sigmas1_n[i][j][k][l] = sigmas1_z[zone_config[i][j]][k][l];
	// 			}
	// 		}			
	// 	}
	// }

    return {Z, x_size, y_size, ref_o, ref_s, ref_l, ref_n, nodes_x, nodes_y, 
    		sigmat_z, sigmas0_z, sigmas1_z, q, h_x, h_y, G, zone_config};

}