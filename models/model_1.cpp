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
void get_model(int &R, int &Z, int &G, int &L, std::vector<double> &h_R, std::vector<std::vector<double>> &q_R, std::vector<std::vector<double>> &sigmat_R, 
      std::vector<std::vector<std::vector<double>>> &sigmas0_R ,std::vector<std::vector<std::vector<double>>> &sigmas1_R, std::vector<std::vector<double>> &sigmat_z, 
      std::vector<std::vector<std::vector<double>>> &sigmas0_z, std::vector<std::vector<std::vector<double>>> &sigmas1_z, 
      std::vector<int> &nodes_R, std::vector<int> &zone_config, int &nxf, double &lb, double &rb){

    R = 1;
    Z = 1;
    G = 1;
    L = 1;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    zone_config = {1};

    sigmat_z = {{1.0}};
    
    sigmas0_z = {{{0.97}}};
    sigmas1_z = {{{0.00}}};

    std::vector<std::vector<double>> q = {{0.0}};

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<double> h_z = {5};

    std::vector<int> nodes_z = {100};

    nxf = 0;
    for (int i = 0 ; i < R ; i++){
        nxf += nodes_z[i];
    }

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    lb = 1;
    rb = 0;

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    h_R.resize(R);
    nodes_R.resize(R);
    sigmat_R.resize(R,vector<double>(G));
    q_R.resize(R,vector<double>(G));
    sigmas0_R.resize(R,vector<vector<double> >(G,vector<double>(G)));
    sigmas1_R.resize(R,vector<vector<double> >(G,vector<double>(G)));

    for(int i = 0 ; i < R ; i++){
        h_R[i] = h_z[i];
        nodes_R[i] = nodes_z[i];
        for(int j = 0 ; j < G ; j++){
            sigmat_R[i][j] = sigmat_z[zone_config[i]-1][j];
            q_R[i][j] = q[i][j];
            for(int k = 0; k < G ; k++){
                sigmas0_R[i][j][k] = sigmas0_z[zone_config[i]-1][j][k];
                sigmas1_R[i][j][k] = sigmas1_z[zone_config[i]-1][j][k];
            }
        }
    }
}