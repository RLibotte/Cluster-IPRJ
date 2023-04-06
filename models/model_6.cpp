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

    const int R = 4;
    const int Z = 4;
    const int G = 4;
    const int L = 1;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<int> zone_config = {1,2,3,4};

    vector<vector <double>> sigmat = {{1.0, 1.0, 1.0, 1.0},{1.0, 1.0, 1.0, 1.0},{1.0, 1.0, 1.0, 1.0},{1.0, 1.0, 1.0, 1.0}};

    // vector<vector<vector<double>>> sigmas0 = {{{0.55, 0.00, 0.00, 0.00},
    //                                            {0.22, 0.60, 0.00, 0.00},
    //                                            {0.10, 0.20, 0.45, 0.00},
    //                                            {0.05, 0.10, 0.25, 0.55}},
                                              
    //                                           {{0.60, 0.00, 0.00, 0.00},
    //                                            {0.15, 0.65, 0.00, 0.00},
    //                                            {0.09, 0.15, 0.70, 0.00},
    //                                            {0.05, 0.07, 0.20, 0.50}},

    //                                            {{0.65, 0.00, 0.00, 0.00},
    //                                             {0.18, 0.75, 0.00, 0.00},
    //                                             {0.08, 0.12, 0.68, 0.00},
    //                                             {0.02, 0.08, 0.22, 0.60}},

    //                                            {{0.45, 0.00, 0.00, 0.00},
    //                                             {0.30, 0.60, 0.00, 0.00},
    //                                             {0.12, 0.18, 0.63, 0.00},
    //                                             {0.07, 0.14, 0.21, 0.60}}};

        vector<vector<vector<double>>> sigmas0 = {{{0.55, 0.22, 0.10, 0.05},
                                                   {0.00, 0.60, 0.20, 0.10},
                                                   {0.00, 0.00, 0.45, 0.25},
                                                   {0.00, 0.00, 0.00, 0.55}},
                                                  
                                                  {{0.60, 0.15, 0.09, 0.05},
                                                   {0.00, 0.65, 0.15, 0.07},
                                                   {0.00, 0.00, 0.70, 0.20},
                                                   {0.00, 0.00, 0.00, 0.50}},

                                                   {{0.65, 0.18, 0.08, 0.02},
                                                    {0.00, 0.75, 0.12, 0.08},
                                                    {0.00, 0.00, 0.68, 0.22},
                                                    {0.00, 0.00, 0.00, 0.60}},

                                                   {{0.45, 0.30, 0.12, 0.07},
                                                    {0.00, 0.60, 0.18, 0.14},
                                                    {0.00, 0.00, 0.63, 0.21},
                                                    {0.00, 0.00, 0.00, 0.60}}};

    // vector<vector<vector<double>>> sigmas0 = {{{0.95, 0.40, 0.10, 0.05},
    //                                            {0.00, 0.55, 0.40, 0.20},
    //                                            {0.00, 0.00, 0.45, 0.25},
    //                                            {0.00, 0.00, 0.00, 0.42}},
                                          
    //                                           {{0.96, 0.30, 0.05, 0.10},
    //                                            {0.00, 0.65, 0.15, 0.18},
    //                                            {0.00, 0.00, 0.70, 0.20},
    //                                            {0.00, 0.00, 0.00, 0.45}},

    //                                           {{0.84, 0.25, 0.09, 0.02},
    //                                            {0.00, 0.75, 0.18, 0.19},
    //                                            {0.00, 0.00, 0.68, 0.22},
    //                                            {0.00, 0.00, 0.00, 0.55}},

    //                                           {{0.98, 0.30, 0.12, 0.07},
    //                                            {0.00, 0.67, 0.23, 0.16},
    //                                            {0.00, 0.00, 0.63, 0.21},
    //                                            {0.00, 0.00, 0.00, 0.51}}};

    vector<vector<vector<double>>> sigmas1 (Z, vector<vector<double>> (G, vector<double> (G, 0)));

    // vector<vector<double>> q = {{0,0,0,0},{1,0.75,0.5,0.25},{0,0,0,0},{0,0,0,0}};
    vector<vector<double>> q = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<double> h = {5,5,5,5};

    int k = 2000;

    vector<int> nodes = {k, k, k, k};

    int nxf = 0;
    for (int i = 0 ; i < R ; i++){
        nxf += nodes[i];
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    double lb = 1; ////////// Use -1 for reflexive boundary condition
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