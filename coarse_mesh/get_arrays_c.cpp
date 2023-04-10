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

tuple<vector<vector<vector<double>>>, vector<vector<double>>> 
				get_arrays_c(int N, int R, int G, vector<vector<vector <double>>> V, vector<vector<double>> ni,
				vector<vector<vector <double>>> VI, vector<vector<double>> ni_i,
				vector<double> hj, vector<vector<vector<double>>> A_p, vector<vector<double>> q){

	const int NG = N*G;
	int r_index = 0;

	vector<double> q_p (NG);
	vector<double> psi_p_temp (NG);

	vector<vector<double>> temp (NG, vector<double> (NG));

	vector<vector<double>> A_p_temp (NG, vector<double> (NG));
	vector<vector<double>> psi_p (R, vector<double> (NG));

	vector<vector<vector<double>>> Am (R, vector<vector<double>> (NG, vector<double> (NG)));

	for(int i = 0 ; i < R ; i++){
		for (int j = 0; j < G; j++) {
	        for (int k = 0; k < N; k++) {
	            if (k < N / 2) {
	                for (int l = 0; l < NG; l++) {
	                    if (ni[l][i] > 0 && ni_i[l][i] == 0) {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l];

	                    } else if (ni[l][i] < 0 && ni_i[l][i] == 0) {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l] * exp(-hj[i] * fabs(ni[l][i]));

	                        ///////////////////////////////////////////////////////////////////////////////////////

	                    } else if (ni[l][i] > 0 && ni_i[l][i] > 0) {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l] * cos(ni_i[l][i] * hj[i] / 2) - VI[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2);

	                        // printf("%.6e \n", cos(ni_i[l][i] * hj[i] / 2) - VI[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2));

	                        printf("%.6e\t--\t%.6e\t--\t%.6e \n", V[i][k + j*N][l], V[i][k + j*N][l] * cos(ni_i[l][i] * hj[i] / 2) - VI[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2), fabs(V[i][k + j*N][l] - V[i][k + j*N][l] * cos(ni_i[l][i] * hj[i] / 2) - VI[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2)));
	                    
	                    } else if (ni[l][i] < 0 && ni_i[l][i] > 0) {
	                        Am[i][k + j*N][l] = exp(-hj[i] * fabs(ni[l][i]))*(V[i][k + j*N][l] * cos(ni_i[l][i]* hj[i] / 2) - VI[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2));

	                        ///////////////////////////////////////////////////////////////////////////////////////

	                    } else if (ni[l][i] > 0 && ni_i[l][i] < 0) {
	                        Am[i][k + j*N][l] = VI[i][k + j*N][l] * cos(ni_i[l][i] * hj[i] / 2) + V[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2);

	                    } else if (ni[l][i] < 0 && ni_i[l][i] < 0) {
	                        Am[i][k + j*N][l] = exp(-hj[i] * fabs(ni[l][i]))*(VI[i][k + j*N][l] * cos(ni_i[l][i] * hj[i] / 2) + V[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2));
	                    }

	                }
	            } else { 
	                for (int l = 0; l < NG; l++) {
	                	if (ni[l][i] > 0 && ni_i[l][i] == 0) {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l] * exp(-hj[i] * fabs(ni[l][i]));

	                    } else if (ni[l][i] < 0 && ni_i[l][i] == 0) {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l];

	                        ///////////////////////////////////////////////////////////////////////////////////////

	                    } else if (ni[l][i] > 0 && ni_i[l][i] > 0) {
	                        Am[i][k + j*N][l] =  exp(-hj[i] * fabs(ni[l][i])) * (V[i][k + j*N][l] * cos(ni_i[l][i] * hj[i] / 2) + VI[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2));
	                    
	                    } else if (ni[l][i] < 0 && ni_i[l][i] > 0) {
	                        Am[i][k + j*N][l] = V[i][k + j*N][l] * cos(ni_i[l][i]* hj[i] / 2) + VI[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2);
	                    
	                        ///////////////////////////////////////////////////////////////////////////////////////
	                    
	                    } else if (ni[l][i] > 0 && ni_i[l][i] < 0) {
	                        Am[i][k + j*N][l] =  exp(-hj[i] * fabs(ni[l][i]))*(VI[i][k + j*N][l] * cos(ni_i[l][i] * hj[i] / 2) - V[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2));

	                    } else if (ni[l][i] < 0 && ni_i[l][i] < 0) {
	                        Am[i][k + j*N][l] = VI[i][k + j*N][l] * cos(ni_i[l][i] * hj[i] / 2) - V[i][k + j*N][l] * sin(ni_i[l][i] * hj[i] / 2);
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

	return {Am, psi_p};
}