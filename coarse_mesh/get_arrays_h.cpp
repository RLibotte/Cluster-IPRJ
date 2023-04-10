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

tuple<vector<vector<vector<double>>>, vector<vector<double>>> 
				get_arrays(int N, int R, int G, vector<vector<vector <double>>> V, vector<vector<double>> ni,
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
	        for (int k = 0; k < N/2; k++) {
                for (int l = 0; l < NG/2; l++) {
                    Am[i][k][l] = V[i][k][l];
                    Am[i][k][l + N/2] = V[i][k + N/2][l];// * exp(-hj[i] / ni[l][i]);

                    Am[i][k + N/2][l] = V[i][k + N/2][l] * exp(-hj[i] / ni[l][i]);
                    Am[i][k + N/2][l + N/2] = V[i][k][l] * exp(-hj[i] / ni[l][i]);
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

    // temp = inverse(A_p_temp, NG);

    // printf("Ap_temp\n");
    // for(int i = 0 ; i < NG ; i++){
    // 	for(int j = 0 ; j < NG ; j++){
    // 		printf("%.6f ", A_p_temp[i][j]);
    // 	}
    // 	printf("\n");
    // }

    // printf("Ap_inv\n");
    // for(int i = 0 ; i < NG ; i++){
    // 	for(int j = 0 ; j < NG ; j++){
    // 		printf("%.6f ", temp[i][j]);
    // 	}
    // 	printf("\n");
    // }

	    psi_p_temp = mult(inverse(A_p_temp, NG), q_p);
	    psi_p[i] = psi_p_temp;

	}


	    printf("Am\n");
	    print3dArray(Am);


	return {Am, psi_p};
}