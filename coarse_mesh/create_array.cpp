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

void  create_arrays(int N, int G, int R, int L, double lb, double rb, std::vector<std::vector<double>> &psi, 
    std::vector<std::vector<double>> &psi_m, std::vector<std::vector<double>> &phi, std::vector<std::vector<double>> &phi_old,
    std::vector<std::vector<double>> &SS, std::vector<std::vector<double>> &alpha, std::vector<double> &alpha_temp, 
    std::vector<double> &avg, std::vector<double> &D_eig, std::vector<double> &DI_eig, std::vector<std::vector<double>> &V_eig,
    std::vector<std::vector<double>> &VI_eig, std::vector<std::vector<double>> &ni, std::vector<std::vector<double>> &ni_i,
    std::vector<std::vector<std::vector <double>>> &V, std::vector<std::vector<std::vector <double>>> &VI){

	const int NG = N*G;

    psi.resize(NG, std::vector<double> (R+1, 0));
    psi_m.resize(R, std::vector<double> (NG, 0));
    SS.resize(G, std::vector<double> (L, 0));
    phi.resize(G, std::vector<double> (R+1, 0));
    phi_old.resize(G, std::vector<double> (R+1, 1));
    alpha.resize(R, std::vector<double> (NG));
    alpha_temp.resize(NG); 
    avg.resize(iter);

    V_eig.resize(NG , std::vector<double> (NG)); 
    VI_eig.resize(NG , std::vector<double> (NG)); 
    D_eig.resize(NG);
    DI_eig.resize(NG);

    ni.resize(NG, std::vector<double> (R));
    ni_i.resize(NG, std::vector<double> (R));

    V.resize(R, std::vector<std::vector<double>> (NG, std::vector<double> (NG)));
    VI.resize(R, std::vector<std::vector<double>> (NG, std::vector<double> (NG)));

    for(int i = 0 ; i < G ; i++){
        for(int j = 0 ; j < N/2 ; j++){
           psi[j][0] = lb; 
           psi[j + N/2][R] = rb;
        }
    }
}