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

#include "../eig/Alglib/linalg.h"
#include "../eig/Alglib/solvers.h"
#include "../eig/Alglib/ap.h"

#include "../solve.cpp"
#include "../print.cpp"

void get_arrays(int N, int R, int G, std::vector<std::vector<std::vector <double>>> &V, std::vector<std::vector<double>> &ni,
	 std::vector<double> &hj, std::vector<std::vector<std::vector<double>>> &A_p, std::vector<std::vector<double>> &q,
	 std::vector<double> &q_p, std::vector<double> &psi_p_temp, std::vector<std::vector<double>> &temp, std::vector<std::vector<double>> &A_p_temp,
	 std::vector<std::vector<double>> &psi_p, std::vector<std::vector<std::vector <double>>> &Am){

	const int NG = N*G;
	int r_index = 0;

	q_p.resize(NG);
	psi_p_temp.resize(NG);
	temp.resize(NG, std::vector<double> (NG));
	A_p_temp.resize(NG, std::vector<double> (NG));
	psi_p.resize(R, std::vector<double> (NG));
	Am.resize(R, std::vector<std::vector<double>> (NG, std::vector<double> (NG)));

	for(int i = 0 ; i < R ; i++){
		for (int j = 0; j < G; j++) {
	        for (int k = 0; k < N; k++) {
	            if (k < N / 2) {
	                for (int l = 0; l < NG; l++) {
	                    if (ni[l][i] > 0) {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l];
	                    } else {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l] * exp(hj[i] / ni[l][i]);
	                    }
	                }
	            } else {
	                for (int l = 0; l < NG; l++) {
	                    if (ni[l][i] > 0) {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l] * exp(-hj[i] / ni[l][i]);
	                    } else {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l];
	                    }
	                }
	            }
	        }
	    }

	    for(int j = 0 ; j < NG ; j++){
	    	for(int k = 0 ; k < NG ; k++){
	    		A_p_temp[j][k] = A_p[j][k][i];
	    	}
	    }
	    
	    for(int j = 0 ; j < G ; j++){
	    	for(int k = 0 ; k < N ; k++){
	    		q_p[k + j*N] = q[i][j];
	    	}
	    }

	    psi_p_temp = mult(inverse(A_p_temp, NG), q_p);
	    psi_p[i] = psi_p_temp;

	}
}