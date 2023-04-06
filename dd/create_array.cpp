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

void create_arrays(int N, int G, int nxf, int L, double lb, double rb, std::vector<std::vector<double>> &psi, 
    std::vector<std::vector<double>> &psi_m, std::vector<std::vector<double>> &phi, std::vector<std::vector<double>> &phi_old,
    std::vector<std::vector<std::vector<double>>> &SS){

	const int NG = N*G;

    psi.resize(NG,std::vector<double>(nxf+1));
    psi_m.resize(NG, std::vector<double> (nxf));
    SS.resize(G, std::vector<std::vector<double>> (nxf, std::vector<double> (2)));
    phi.resize(G, std::vector<double> (nxf + 1, 0));
    phi_old.resize(G, std::vector<double> (nxf + 1, 1));

    for(int i = 0 ; i < G ; i++){
        for(int j = 0 ; j < N/2 ; j++){
           psi[j][0] = lb;
           psi[j + N/2][nxf] = rb;
        }
    }
}