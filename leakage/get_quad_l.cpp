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

void quad_l_x(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_y, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, int lb, int ub){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// OESTE

	if (index_x == 0){

		if(index_y == 0){ ///////////////////////////////// NOROESTE

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_x[index_y][index_x][aux] + psi_x[index_y + 1][index_x][aux] + psi_y[index_y + 1][index_x][aux])/3;
					p2 = (psi_x[index_y][index_x][aux] + psi_y[index_y][index_x][aux])/2;

					p2 = (nhy[index_y][index_x][k]) * (p2 - p1);
					p1 = (h_x[index_y][index_x + 1] * l_y[index_y][index_x][aux] + h_x[index_y][index_x] * l_y[index_y][index_x + 1][aux])/(h_x[index_y][index_x] + h_x[index_y][index_x + 1]);

					c1[index_y][index_x][aux] = 0.5 * (p1 - p2);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_y[index_y][index_x][aux]);
				}
			}

		} else if (index_y == (nyf - 1)){ ///////////////////////////////// SUDOESTE

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_x[index_y][index_x][aux] + psi_y[index_y + 1][index_x][aux])/2;
					p2 = (psi_x[index_y - 1][index_x][aux] + psi_x[index_y][index_x][aux] + psi_y[index_y][index_x][aux])/3;

					p2 = (nhy[index_y][index_x][k]) * (p2 - p1);
					p1 = (h_x[index_y][index_x + 1] * l_y[index_y][index_x][aux] + h_x[index_y][index_x] * l_y[index_y][index_x + 1][aux])/(h_x[index_y][index_x + 1] + h_x[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p1 - p2);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_y[index_y][index_x][aux]);
				}
			}

		} else { ////////////////////////////// MEIO

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_x[index_y][index_x][aux] + psi_x[index_y + 1][index_x][aux] + psi_y[index_y + 1][index_x][aux])/3;
					p2 = (psi_x[index_y - 1][index_x][aux] + psi_x[index_y][index_x][aux] + psi_y[index_y][index_x][aux])/3;

					p2 = (nhy[index_y][index_x][k]) * (p2 - p1);
					p1 = (h_x[index_y][index_x+1] * l_y[index_y][index_x][aux] + h_x[index_y][index_x] * l_y[index_y][index_x + 1][aux])/(h_x[index_y][index_x+1] + h_x[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p1 - p2);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_y[index_y][index_x][aux]);
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// LESTE

	 else if (index_x == nxf - 1){

		if(index_y == 0){ ///////////////////////////////// NORDESTE

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_x[index_y][index_x + 1][aux] + psi_x[index_y + 1][index_x + 1][aux] + psi_y[index_y + 1][index_x][aux])/3;
					p2 = (psi_x[index_y][index_x + 1][aux] + psi_y[index_y][index_x][aux])/2;

					p2 = (nhy[index_y][index_x][k]) * (p2 - p1);
					p1 = (h_x[index_y][index_x] * l_y[index_y][index_x - 1][aux] + h_x[index_y][index_x - 1] * l_y[index_y][index_x][aux])/(h_x[index_y][index_x - 1] + h_x[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_y[index_y][index_x][aux]);
				}
			}

		} else if (index_y == (nyf - 1)){ ///////////////////////////////// SUDESTE

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_x[index_y][index_x + 1][aux] + psi_y[index_y + 1][index_x][aux])/2;
					p2 = (psi_x[index_y - 1][index_x + 1][aux] + psi_x[index_y][index_x + 1][aux] + psi_y[index_y][index_x][aux])/3;

					p2 = (nhy[index_y][index_x][k]) * (p2 - p1);
					p1 = (h_x[index_y][index_x] * l_y[index_y][index_x - 1][aux] + h_x[index_y][index_x - 1] * l_y[index_y][index_x][aux])/(h_x[index_y][index_x - 1] + h_x[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_y[index_y][index_x][aux]);
				}
			}

		} else { ////////////////////////////// MEIO

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_x[index_y][index_x + 1][aux] + psi_x[index_y + 1][index_x + 1][aux] + psi_y[index_y + 1][index_x][aux])/3;
					p2 = (psi_x[index_y][index_x + 1][aux] + psi_x[index_y - 1][index_x + 1][aux] + psi_y[index_y][index_x][aux])/3;

					p2 = (nhy[index_y][index_x][k]) * (p2 - p1);
					p1 = (h_x[index_y][index_x - 1] * l_y[index_y][index_x][aux] + h_x[index_y][index_x] * l_y[index_y][index_x - 1][aux])/(h_x[index_y][index_x - 1] + h_x[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_y[index_y][index_x][aux]);
				}
			}
		}

	} else {

		for(int g = 0 ; g < G ; g++){
			for(int k = lb ; k < ub ; k++){
				aux = k + g*M;

				p1 = (h_x[index_y][index_x + 1] * l_y[index_y][index_x][aux] + h_x[index_y][index_x] * l_y[index_y][index_x + 1][aux])/(h_x[index_y][index_x + 1] + h_x[index_y][index_x]);
				p2 = (h_x[index_y][index_x - 1] * l_y[index_y][index_x][aux] + h_x[index_y][index_x] * l_y[index_y][index_x - 1][aux])/(h_x[index_y][index_x - 1] + h_x[index_y][index_x]);

				c1[index_y][index_x][aux] = 0.5 * (p1 - p2);
				c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_y[index_y][index_x][aux]);
			}
		}
	}

	// for(int g = 0 ; g < G ; g++){
	// 	for(int k = lb ; k < ub ; k++){
	// 		aux = k + g*M;
	// 		c1[index_y][index_x][aux] = 0.5 * (p1 - p2);
	// 		c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_y[index_y][index_x][aux]);
	// 	}
	// }

	b2[index_y][index_x] = mult(A_p_inv_z[zcn[index_y][index_x]], c2[index_y][index_x]);

	for(int k = 0 ; k < M ; k++){
		for(int g = 0 ; g < G ; g++){
			aux = k + g*M;

			arr[aux] = 6 * mhx[index_y][index_x][k] * b2[index_y][index_x][aux] - c1[index_y][index_x][aux];
		}
	}

	b1[index_y][index_x] = mult(A_p_inv_z[zcn[index_y][index_x]], arr);

	for(int k = 0 ; k < M ; k++){
		for(int g = 0 ; g < G ; g++){
			aux = k + g*M;
			arr[aux] = q[index_y][index_x][g] - l_y[index_y][index_x][aux] - 2 *(mhx[index_y][index_x][k]) * b1[index_y][index_x][aux];
		}
	}

	b0[index_y][index_x] = mult(A_p_inv_z[zcn[index_y][index_x]], arr);

	// for(int i = 0 ; i < M*G ; i++){
	// 	b_aux_x[index_y][index_x][i] = b0[index_y][index_x][i] - b2[index_y][index_x][i];
	// }

}


// -------------------------------------------------------------------------------------------------------------------------------------------------------------


void quad_l_y(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_x, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, int lb, int ub){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// NORTE

	if(index_y == 0) {

		if(index_x == 0) { ////////////////////////////////////////////////// NOROESTE

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_y[index_y][index_x][aux] + psi_x[index_y][index_x][aux])/2;
					p2 = (psi_y[index_y][index_x][aux] + psi_y[index_y][index_x+1][aux] + psi_x[index_y][index_x+1][aux])/3;

					p2 = (mhx[index_y][index_x][k]) * (p2 - p1);
					p1 = (h_y[index_y+1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y + 1][index_x][aux])/(h_y[index_y+1][index_x] + h_y[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_x[index_y][index_x][aux]);
				}
			}

		} else if(index_x == nxf-1) { ////////////////////////////////////////////////// NORDESTE

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_y[index_y][index_x][aux] + psi_x[index_y][index_x+1][aux])/2;
					p2 = (psi_y[index_y][index_x][aux] + psi_y[index_y][index_x-1][aux] + psi_x[index_y][index_x][aux])/3;

					p2 = (mhx[index_y][index_x][k]) * (p1 - p2);
					p1 = (h_y[index_y+1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y + 1][index_x][aux])/(h_y[index_y][index_x] + h_y[index_y+1][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_x[index_y][index_x][aux]);
				}
			}

		} else { ////////////////////////////////////////////////// MEIO

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_y[index_y][index_x-1][aux] + psi_y[index_y][index_x][aux] + psi_x[index_y][index_x][aux])/3;
					p2 = (psi_y[index_y][index_x][aux] + psi_y[index_y][index_x+1][aux] + psi_x[index_y][index_x+1][aux])/3;

					p2 = (mhx[index_y][index_x][k]) * (p2 - p1);
					p1 = (h_y[index_y+1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y + 1][index_x][aux])/(h_y[index_y+1][index_x] + h_y[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_x[index_y][index_x][aux]);
				}
			}
		}
	
	} else if (index_y == nyf - 1){ //////////////////////////////////////////////////////////////////////////////////////////////////// SUL

		if(index_x == 0){ ////////////////////////////////////////////////// SUDOESTE

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_y[index_y + 1][index_x][aux] + psi_x[index_y][index_x][aux])/2;
					p2 = (psi_y[index_y + 1][index_x][aux] + psi_y[index_y+1][index_x + 1][aux] + psi_x[index_y][index_x + 1][aux])/3;

					p1 = (mhx[index_y][index_x][k]) * (p2 - p1);
					p2 = (h_y[index_y-1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y - 1][index_x][aux])/(h_y[index_y-1][index_x] + h_y[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_x[index_y][index_x][aux]);
				}
			}
		
		} else if(index_x == nxf - 1){ ////////////////////////////////////////////////// SUDESTE

		    for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_y[index_y + 1][index_x][aux] + psi_x[index_y][index_x + 1][aux])/2;
					p2 = (psi_y[index_y + 1][index_x - 1][aux] + psi_y[index_y + 1][index_x][aux] + psi_x[index_y][index_x][aux])/3;

					p1 = (mhx[index_y][index_x][k]) * (p1 - p2);
					p2 = (h_y[index_y - 1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y - 1][index_x][aux])/(h_y[index_y - 1][index_x] + h_y[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_x[index_y][index_x][aux]);
				}
			}

		} else { ////////////////////////////////////////////////// MEIO

			for(int g = 0 ; g < G ; g++){
				for(int k = lb ; k < ub ; k++){
					aux = k + g*M;

					p1 = (psi_y[index_y + 1][index_x - 1][aux] + psi_y[index_y + 1][index_x][aux] + psi_x[index_y][index_x][aux])/3;
					p2 = (psi_y[index_y + 1][index_x][aux] + psi_y[index_y + 1][index_x + 1][aux] + psi_x[index_y][index_x + 1][aux])/3;

					p1 = (mhx[index_y][index_x][k]) * (p2 - p1);
					p2 = (h_y[index_y - 1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y - 1][index_x][aux])/(h_y[index_y - 1][index_x] + h_y[index_y][index_x]);

					c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
					c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_x[index_y][index_x][aux]);
				}
			}
		}

	} else {

	 	for(int g = 0 ; g < G ; g++){
			for(int k = lb ; k < ub ; k++){
				aux = k + g*M;

				p1 = (h_y[index_y+1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y + 1][index_x][aux])/(h_y[index_y+1][index_x] + h_y[index_y][index_x]);
				p2 = (h_y[index_y-1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y - 1][index_x][aux])/(h_y[index_y-1][index_x] + h_y[index_y][index_x]);

				c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
				c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_x[index_y][index_x][aux]);
			}
		}
	}

	b2[index_y][index_x] = mult(A_p_inv_z[zcn[index_y][index_x]], c2[index_y][index_x]);

	for(int k = 0 ; k < M ; k++){
		for(int g = 0 ; g < G ; g++){
			aux = k + g*M;
			arr[aux] = 6 * nhy[index_y][index_x][k] * b2[index_y][index_x][aux] - c1[index_y][index_x][aux];
		}
	}

	b1[index_y][index_x] = mult(A_p_inv_z[zcn[index_y][index_x]], arr);

	for(int k = 0 ; k < M ; k++){
		for(int g = 0 ; g < G ; g++){
			aux = k + g*M;
			arr[aux] = q[index_y][index_x][g] - l_x[index_y][index_x][aux] - 2 * nhy[index_y][index_x][k] * b1[index_y][index_x][aux];
		}
	}

	b0[index_y][index_x] = mult(A_p_inv_z[zcn[index_y][index_x]], arr);

	// for(int i = 0 ; i < M*G ; i++){
	// 	b_aux_y[index_y][index_x][i] = b0[index_y][index_x][i] - b2[index_y][index_x][i];
	// }
}