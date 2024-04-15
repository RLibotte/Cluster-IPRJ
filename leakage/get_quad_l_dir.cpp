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

void get_b(int index_x, int index_y, int M, int G, vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l, 
			vector<vector<vector<double>>> &q, vector<vector<vector<double>>> &c1, vector<vector<vector<double>>> &c2, 
			vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1, vector<vector<vector<double>>> &b2,
			vector<vector<vector<double>>> &mn, vector<vector<vector<double>>> &b_aux_x, vector<double> &arr){

	int aux;

	b2[index_y][index_x] = mult(A_p_inv_z[zcn[index_y][index_x]], c2[index_y][index_x]);

	for(int k = 0 ; k < M ; k++){
		for(int g = 0 ; g < G ; g++){
			aux = k + g*M;

			arr[aux] = 6 * mn[index_y][index_x][k] * b2[index_y][index_x][aux] - c1[index_y][index_x][aux];
		}
	}

	b1[index_y][index_x] = mult(A_p_inv_z[zcn[index_y][index_x]], arr);

	for(int k = 0 ; k < M ; k++){
		for(int g = 0 ; g < G ; g++){
			aux = k + g*M;
			arr[aux] = q[index_y][index_x][g] - l[index_y][index_x][aux] - 2 *(mn[index_y][index_x][k]) * b1[index_y][index_x][aux];
		}
	}

	b0[index_y][index_x] = mult(A_p_inv_z[zcn[index_y][index_x]], arr);

	for(int i = 0 ; i < M*G ; i++){
		b_aux_x[index_y][index_x][i] = b0[index_y][index_x][i] - b2[index_y][index_x][i];
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void quad_x_sw(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_y, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_x){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_y, q, c1, c2, b0, b1, b2, mhx, b_aux_x, arr);

}

void quad_x_nw(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_y, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_x){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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


	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_y, q, c1, c2, b0, b1, b2, mhx, b_aux_x, arr);

}

void quad_x_w(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_y, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_x){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_y, q, c1, c2, b0, b1, b2, mhx, b_aux_x, arr);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void quad_x_se(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_y, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_x){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_y, q, c1, c2, b0, b1, b2, mhx, b_aux_x, arr);

}

void quad_x_ne(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_y, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_x){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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


	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_y, q, c1, c2, b0, b1, b2, mhx, b_aux_x, arr);

}

void quad_x_e(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_y, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_x){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_y, q, c1, c2, b0, b1, b2, mhx, b_aux_x, arr);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void quad_x(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_y, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_x){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

	for(int g = 0 ; g < G ; g++){
		for(int k = lb ; k < ub ; k++){
			aux = k + g*M;

			p1 = (h_x[index_y][index_x + 1] * l_y[index_y][index_x][aux] + h_x[index_y][index_x] * l_y[index_y][index_x + 1][aux])/(h_x[index_y][index_x + 1] + h_x[index_y][index_x]);
			p2 = (h_x[index_y][index_x - 1] * l_y[index_y][index_x][aux] + h_x[index_y][index_x] * l_y[index_y][index_x - 1][aux])/(h_x[index_y][index_x - 1] + h_x[index_y][index_x]);

			c1[index_y][index_x][aux] = 0.5 * (p1 - p2);
			c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_y[index_y][index_x][aux]);
		}
	}

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_y, q, c1, c2, b0, b1, b2, mhx, b_aux_x, arr);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// -------------------------------------------------------------------------------------------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void quad_y_nw(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_x, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_y){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_x, q, c1, c2, b0, b1, b2, nhy, b_aux_y, arr);

}

void quad_y_ne(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_x, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_y){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_x, q, c1, c2, b0, b1, b2, nhy, b_aux_y, arr);

}

void quad_y_n(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_x, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_y){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_x, q, c1, c2, b0, b1, b2, nhy, b_aux_y, arr);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void quad_y_sw(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_x, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_y){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_x, q, c1, c2, b0, b1, b2, nhy, b_aux_y, arr);

}

void quad_y_se(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_x, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_y){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_x, q, c1, c2, b0, b1, b2, nhy, b_aux_y, arr);

}

void quad_y_s(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_x, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_y){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

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

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_x, q, c1, c2, b0, b1, b2, nhy, b_aux_y, arr);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void quad_y(int index_x, int index_y, int M, int G, int nxf, int nyf, vector<double> &mi, vector<double> &n,
			vector<vector<int>> &zcn, vector<vector<vector<double>>> &A_p_inv_z, vector<vector<vector<double>>> &l_x, 
			vector<vector<vector<double>>> &q, vector<vector<double>> &h_x, vector<vector<double>> &h_y, vector<vector<vector<double>>> &c1, 
			vector<vector<vector<double>>> &c2, vector<vector<vector<double>>> &b0, vector<vector<vector<double>>> &b1,
			vector<vector<vector<double>>> &b2, vector<vector<vector<double>>> &psi_x, vector<vector<vector<double>>> &psi_y,
			int lb, int ub, vector<vector<vector<double>>> &mhx, vector<vector<vector<double>>> &nhy, 
			vector<vector<vector<double>>> &b_aux_y){

	const int MG = M*G;

	int aux;
	double p1, p2;
	vector<double> arr (MG);

	for(int g = 0 ; g < G ; g++){
		for(int k = lb ; k < ub ; k++){
			aux = k + g*M;

			p1 = (h_y[index_y+1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y + 1][index_x][aux])/(h_y[index_y+1][index_x] + h_y[index_y][index_x]);
			p2 = (h_y[index_y-1][index_x] * l_x[index_y][index_x][aux] + h_y[index_y][index_x] * l_x[index_y - 1][index_x][aux])/(h_y[index_y-1][index_x] + h_y[index_y][index_x]);

			c1[index_y][index_x][aux] = 0.5 * (p2 - p1);
			c2[index_y][index_x][aux] = -(0.5 * (p1 + p2) - l_x[index_y][index_x][aux]);
		}
	}

	get_b(index_x, index_y, M, G, zcn, A_p_inv_z, l_x, q, c1, c2, b0, b1, b2, nhy, b_aux_y, arr);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////