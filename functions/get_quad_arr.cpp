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

 tuple <vector<int>,
 		vector<int>,
 		vector<int>,
 		vector<int>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>,
 		vector<vector<vector<double>>>> 

 		get_quad_arr(int nxf, int nyf, int MG){

 	vector<vector<vector<double>>> c1_x (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   c2_x (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   b0_x (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   b1_x (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   b2_x (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   b02_x (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   
								   c1_y (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   c2_y (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   b0_y (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   b1_y (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   b2_y (nyf, vector<vector<double>> (nxf, vector <double>(MG))),
								   b02_y (nyf, vector<vector<double>> (nxf, vector <double>(MG)));


		vector<int> nxf_arr(nxf), 
					nxf_arr_inv(nxf), 
					nyf_arr(nyf),
					nyf_arr_inv(nyf); 

	   for(int i = 0 ; i < nxf ; i++){
	   		nxf_arr[i] = i;
	   		nxf_arr_inv[i] = nxf - i - 1;
	    }

	   	for(int i = 0 ; i < nyf ; i++){
	   		nyf_arr[i] = i;
	   		nyf_arr_inv[i] = nyf - i - 1;
	    }

	    return {nxf_arr, nyf_arr, nxf_arr_inv, nyf_arr_inv, 
	    c1_x, c2_x, b0_x, b1_x, b2_x, b02_x, 
	    c1_y, c2_y, b0_y, b1_y, b2_y, b02_y};
 }