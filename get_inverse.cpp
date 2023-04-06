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

#include "./eig/Alglib/linalg.h"
#include "./eig/Alglib/solvers.h"
#include "./eig/Alglib/ap.h"

tuple <vector<vector<vector<double>>>> get_inverse(const int NG, int R, vector<vector<vector<double>>> arr){

	vector<vector<double>> arr_temp (NG, vector<double> (NG));
	vector<vector<vector<double>>> ans (R, vector<vector<double>> (NG, vector<double> (NG)));

	for(int i = 0 ; i < R ; i++){
		ans[i] = inverse(arr[i], NG);
	}

	return {ans};

}

