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

std::vector<std::vector<double>> sweep_lr(std::vector<std::vector<double>> psi, std::vector<double> mi, std::vector<std::vector<double>> sigmat,
    std::vector<std::vector<std::vector<double>>> SS, std::vector<std::vector<double>> q, std::vector<double> h, int nxf, int G, int N){
	for(int k = 0 ; k < nxf ; k++){    
        // #pragma omp parallel for collapse(2)
        for(int g = 0 ; g < G ; g++){
            for(int j = 0 ; j < N/2 ; j++){
                psi[j + g*N][k+1] = (psi[j + g*N][k] * (mi[j]/h[k] - sigmat[k][g]/2) + SS[g][k][0] 
                    + 1.5 * SS[g][k][1] * mi[j] + q[k][g]) / (mi[j]/h[k] + sigmat[k][g]/2);
            }
        }
    }
    return psi;
}

std::vector<std::vector<double>> sweep_rl(std::vector<std::vector<double>> psi, std::vector<double> mi, std::vector<std::vector<double>> sigmat,
    std::vector<std::vector<std::vector<double>>> SS, std::vector<std::vector<double>> q, std::vector<double> h, int nxf, int G, int N){
	for(int k = nxf-1 ; k >= 0 ; k--){
        // #pragma omp parallel for collapse(2)
        for(int g = 0 ; g < G ; g++){
            for(int j = N/2 ; j < N ; j++){
                psi[j + g*N][k] =  (psi[j + g*N][k+1] * (fabs(mi[j])/h[k] - sigmat[k][g]/2) + SS[g][k][0] 
                    + 1.5 * SS[g][k][1] * mi[j] + q[k][g])/ (fabs(mi[j])/h[k] + sigmat[k][g]/2);
            }
        }
    }
    return psi;
}