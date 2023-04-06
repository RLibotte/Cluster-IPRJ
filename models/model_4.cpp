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

tuple<int, int, int, int, vector<double>, vector<vector<double>>, vector<vector<double>>, 
      vector<vector<vector<double>>>,vector<vector<vector<double>>>, vector<vector<double>>, 
      vector<vector<vector<double>>>,vector<vector<vector<double>>>, 
      vector<int>, vector<int>, int, double, double> get_model(){

    const int R = 1;
    const int Z = 1;
    const int G = 2;
    const int L = 1;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<int> zone_config = {1};

    vector<vector <double>> sigmat = {{1.0, 1.0}};

    vector<vector<vector<double>>> sigmas0 = {{{0.99, 0.008},
                                               {0.005,0.98}}};
                                               
    vector<vector<vector<double>>> sigmas1 (Z, vector<vector<double>> (G, vector<double> (G, 0)));

    vector<vector<double>> q = {{0}};

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<double> h = {100};

    vector<int> nodes = {10000};

    int nxf = 0;
    for (int i = 0 ; i < R ; i++){
        nxf += nodes[i];
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    double lb = 1;
    double rb = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<double> h_R (R);
    vector<int> nodes_R (R);
    vector<vector <double>> sigmat_R(R, vector<double> (G));
    vector<vector<vector<double>>> sigmas0_R (R, vector<vector<double>> (G, vector<double> (G)));
    vector<vector<vector<double>>> sigmas1_R (R, vector<vector<double>> (G, vector<double> (G)));
    vector<vector<double>> q_R (R, vector<double> (G));

    for(int i = 0 ; i < R ; i++){
        h_R[i] = h[i];
        nodes_R[i] = nodes[i];
        for(int j = 0 ; j < G ; j++){
            sigmat_R[i][j] = sigmat[zone_config[i]-1][j];
            q_R[i][j] = q[i][j];
            for(int k = 0; k < G ; k++){
                sigmas0_R[i][j][k] = sigmas0[zone_config[i]-1][j][k];
                sigmas1_R[i][j][k] = sigmas1[zone_config[i]-1][j][k];
            }
        }
    }
    
    return {R, Z, G, L, h_R, q_R, sigmat_R, sigmas0_R, sigmas1_R, sigmat, sigmas0, sigmas1, nodes_R, zone_config, nxf, lb, rb};
}