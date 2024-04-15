#include <iostream>
#include <vector>

tuple<double, double, double, double> getS(const int M, const int G, vector<double> &w, vector<double> &mi, vector<double> &n,
		  					 vector<double> psi_x_b, vector <double> psi_y_b){

	double S0_x = 0, S0_y = 0, S1_x = 0, S1_y = 0;

	// for(int i = 0 ; i < G ; i++){
 //        for(int j = 0 ; j < G ; j++){
 //            for(int l = 0 ; l < N ; l++){
 //                SS[i] = SS[i] + (sigmas0[j][i][k]/2)*(w[l] * psi_bar[r_index][k]);
 //                r_index++;
 //            }
 //        }
 //        r_index = 0;
 //    }

	for(int k = 0 ; k < M ; k++){
		S0_x = S0_x + w[k] * psi_x_b[k];
		S0_y = S0_y + w[k] * psi_y_b[k];

		S1_x = S1_x + mi[k] * w[k] * psi_x_b[k];
		S1_y = S1_y + n[k] * w[k] * psi_y_b[k];	
	}

	S1_x = S1_x; 
	S1_y = S1_y;

	return {S0_x, S0_y, S1_x, S1_y};
}