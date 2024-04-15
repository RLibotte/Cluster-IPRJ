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

tuple<vector<vector<vector<double>>>, 
	  vector<vector<vector<double>>>,
	  vector<vector<vector<double>>>,
	  vector<vector<vector<double>>>,
	  vector<vector<vector<double>>>,
	  vector<vector<vector<double>>>,
	  vector<vector<vector<double>>>,
	  vector<vector<vector<vector<double>>>>,
	  vector<vector<vector<vector<double>>>>,
	  vector<vector<vector<double>>>,
	  vector<vector<vector<double>>>,
	  vector<vector<vector<double>>>>
  	
  	get_matrix(int x_size, int y_size, int nxf, int nyf, int M, int G){

  	const int MG = M*G;

	vector<vector<vector<double>>> psi_x_b (nyf, vector<vector<double> > (nxf, vector <double>(MG))),
								   psi_y_b (nyf, vector<vector<double> > (nxf, vector <double>(MG))),
								   scalar_flux_x (nyf, vector<vector<double> > (nxf, vector <double>(G))),
							 	   scalar_flux_x_old (nyf, vector<vector<double> > (nxf, vector <double>(G))),
								   scalar_flux_x_out (y_size, vector<vector<double> > (x_size, vector <double>(G))),
								   SS_x (nyf, vector<vector<double> > (nxf, vector <double>(G))),
   								   SS_y (nyf, vector<vector<double> > (nxf, vector <double>(G))),
   								   psi_x (nyf, vector<vector<double> > (nxf + 1, vector <double>(MG))),
   								   psi_x_old (nyf, vector<vector<double> > (nxf + 1, vector <double>(MG))),
   								   psi_y (nyf + 1, vector<vector<double> > (nxf, vector <double>(MG)));


	vector<vector<vector<vector<double>>>> SS_x_1 (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(G, vector<double> (M)))),
										   SS_y_1 (nyf, vector<vector<vector<double>>> (nxf, vector<vector<double>>(G, vector<double> (M))));

	return {psi_x_b, psi_y_b,
    		scalar_flux_x, scalar_flux_x_old, scalar_flux_x_out, 
    		SS_x, SS_y, SS_x_1, SS_y_1, psi_x, psi_x_old, psi_y};
}