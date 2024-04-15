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

tuple<vector<double>, vector<double>, vector<double>, vector<int>> getLQn(int N){

	int M = N*(N + 2)/2;

    vector<double> mi (M);
    vector<double> n (M);
    vector<double> w (M);
    vector<int> order (M);
    
	string filename = "../../LQn/S" + to_string(N) + ".txt";

	vector<double> data(4 * M); 
	std::ifstream ifile(filename, std::ios::in);
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input Quadrature file!\n";
    }
    double num = 0.0;
    int k = 0 ;
    while (ifile >> num) {
        data[k] = num;
        k++;
    }

    for(int i = 0 ; i < M ; i++){
        mi[i] = data[i];
        n[i] = data[i + M];
        w[i] = data[i + 2 * M];
        order[i] = data[i + 3 * M];
    }

	return {mi, n, w, order};
}