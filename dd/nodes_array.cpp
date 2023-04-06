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



void nodes_arr(int N, int L, int R, int G, int nxf, std::vector<int> &nodes, 
  std::vector<std::vector<double>> &sigmat_R, std::vector<std::vector<std::vector <double>>> &sigmas0_R,
  std::vector<std::vector<std::vector <double>>> &sigmas1_R, std::vector<std::vector <double>> &q_R, 
  std::vector<double> &h_R,  std::vector<std::vector<double>> &sigmat, std::vector<std::vector<double>> &q, 
  std::vector<std::vector<std::vector<double>>> &sigmas0, std::vector<std::vector<std::vector<double>>> &sigmas1, 
  std::vector<double> &h, std::vector<std::vector <double>> &psi_m, std::vector<std::vector<std::vector <double>>> &SS,
  std::vector<double> &pos_arr){

    h.resize(nxf);
    pos_arr.resize(nxf+1);

    sigmat.resize(nxf, std::vector<double> (G));
    q.resize(nxf, std::vector<double> (G));
    psi_m.resize(N*G, std::vector<double> (nxf));

    sigmas0.resize(nxf, std::vector<std::vector<double>> (G, std::vector<double> (G)));
    sigmas1.resize(nxf, std::vector<std::vector<double>> (G, std::vector<double> (G)));
    SS.resize(G, std::vector<std::vector<double>> (nxf, std::vector<double> (L)));

    int index = 0;
    for(int i = 0 ; i < R ; i++){
        for(int k = 0 ; k < nodes[i] ; k++){
            h[index] = h_R[i]/nodes[i];
            for(int j = 0 ; j < G ; j++){
                sigmat[index][j] = sigmat_R[i][j];
                q[index][j] = q_R[i][j];
                for (int l = 0 ; l < G ; l++){
                    sigmas0[index][j][l] = sigmas0_R[i][j][l];
                    sigmas1[index][j][l] = sigmas1_R[i][j][l];
                }
            }
            index++;
        }
    }
}