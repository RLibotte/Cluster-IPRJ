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

    const int R = 3;
    const int Z = 2;
    const int G = 1;
    const int L = 2;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<int> zone_config = {1, 2, 1};

    vector<vector <double>> sigmat = {{1.0}, {0.6}};
    
    vector<vector<vector<double>>> sigmas0 = {{{0.9}},{{0.4}}};
    vector<vector<vector<double>>> sigmas1 = {{{0.8}},{{0.3}}};
    
    vector<vector<double>> q = {{0},{0},{0}};

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<double> h = {20, 50, 30};

    int cte = 100;

    vector<int> nodes = {20 * cte, 50 * cte, 30 * cte};

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