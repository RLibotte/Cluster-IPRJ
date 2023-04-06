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

tuple<vector<double>, vector<vector<double>>> 
	print_ref(int N, int G, int NG, int R, int nxf, double lb, double rb, 
		vector<double> h, vector<int> nodes, vector<vector<double>> alpha, 
		vector<vector<vector<double>>> V, vector<vector<double>> ni, vector<vector<double>> psi_p,
		vector<double> w, vector<vector<double>> phi){

	double total = 0;
    double pos = 0;

    vector<vector<double>> psi_nxf(NG, vector<double> (nxf + 1, 0));
    vector<vector<double>> phi_nxf(G, vector<double> (nxf + 1, 0)); 
    vector<double> pos_arr (nxf+1);
    vector<double> pos_r (nxf+1);

    int ind = 1;
    for(int i = 0 ; i < R ; i++){
        for(int j = 0 ; j < nodes[i] ; j++){
            pos_arr[ind] = pos_arr[ind - 1] + h[i]/nodes[i];
            ind++;
        }
    }

    double pos_temp = 0;
    ind = 1;
    for(int i = 0 ; i < R; i++){
        for(int j = 0 ; j < nodes[i] - 1; j++){
            pos_r[ind] = pos_temp + h[i]/nodes[i];
            pos_temp += h[i]/nodes[i];
            ind++;
        }
        pos_temp = 0;
        ind++;
    }
    pos_r[ind - 1] = pos_r[ind - 2] + h[R-1]/nodes[R-1];

    for(int i = 0 ; i < G ; i++){
        for(int j = 0 ; j < N/2 ; j++){
           psi_nxf[j][0] = lb;
           psi_nxf[j + N/2][nxf] = rb;
        }
    }

    ind = 0;

    for(int i = 0 ; i < R ; i++){
        for(int j = 0 ; j < nodes[i]; j++){
            for(int g = 0 ; g < G ; g++){
                for(int k = 0 ; k < N ; k++){
                    if (k < N / 2) {
                        for (int l = 0; l < NG; l++) {
                            if (ni[l][i] > 0) { 
                                psi_nxf[k + g*N][ind+1] += alpha[i][l] * V[i][k + g*N][l] * exp(-pos_r[ind+1]/ni[l][i]) + psi_p[i][k + g*N]/NG;
                            } else {
                                psi_nxf[k + g*N][ind+1] += alpha[i][l] * V[i][k + g*N][l] * exp(-(pos_r[ind+1] - h[i])/ni[l][i]) + psi_p[i][ + g*N]/NG;
                            }
                        }
                    } else {
                        for (int l = 0; l < N * G; l++) {
                            if (ni[l][i] > 0) {
                                psi_nxf[k + g*N][ind] += alpha[i][l] * V[i][k + g*N][l] * exp(-(pos_r[ind])/ni[l][i]) + psi_p[i][k + g*N]/NG;
                            } else {
                                psi_nxf[k + g*N][ind] += alpha[i][l] * V[i][k + g*N][l] * exp(-(pos_r[ind] - h[i])/ni[l][i]) + psi_p[i][k + g*N]/NG;
                            }
                        }
                    }
                }
            }
            ind++;
        }
    }

    for(int i = 0 ; i < nxf + 1; i++){
        for(int g = 0 ; g < G ; g++){
            for(int k = 0 ; k < N ; k++){
                phi_nxf[g][i] += w[k] * psi_nxf[k + g*N][i];
            }
            phi_nxf[g][i] = phi_nxf[g][i]/2; 
        }
    }

    ind = 0;
    for(int i = 0 ; i < R ; i++){
        ind += nodes[i];
        for(int g = 0 ; g < G ; g++){
            phi_nxf[g][ind] = phi[g][i+1];
        }
    }
    for(int g = 0 ; g < G ; g++){
        phi_nxf[g][0] = phi[g][0];
    }

    return {pos_arr, phi_nxf};
}