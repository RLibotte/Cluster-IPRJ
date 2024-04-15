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

void get_psi_out_3(int j,
			 int i, 
			 const int &M, 
			 const int &G, 
			 int &nxf, 
			 int &nyf, 
			 vector<vector<vector<double>>> &l_x,
	 		 vector<vector<vector<double>>> &l_y, 
	 		 vector<double> &mi, 
	 		 vector<double> &n, 
	 		 vector<double> &w,
			 vector<vector<vector<double>>> &psi_x, 
			 vector<vector<vector<double>>> &psi_y, 
			 vector<vector<vector<double>>> &q,
			 vector<double> &b_in_x, 
			 vector<double> &b_in_y, 
			 vector<vector<double>> &h_x, 
			 vector<vector<double>> &h_y,
			 vector<double> &psi_x_out, 
			 vector<double> &psi_y_out, 
			 vector<vector<vector<double>>> &b1_x, 
			 vector<vector<vector<double>>> &b1_y,
			 vector<vector<int>> &zcn, 
			 vector<vector<vector<double>>> &mhx, 
			 vector<vector<vector<double>>> &nhy, 
			 vector<vector<vector<double>>> &b_aux_x, 
			 vector<vector<vector<double>>> &b_aux_y,
			 vector<vector<vector<vector<double>>>> &psi_x_aux,
			 vector<vector<vector<vector<double>>>> &psi_y_aux){

	int index;

	get_b_in(i, j, G, M, b_in_x, b_in_y, b_aux_x, b1_x, b_aux_y, b1_y, psi_x, psi_y);
			
	mult_part(psi_x_aux[i][j], b_in_x, psi_x_out, M/2, 3*M/4, G, M);
	mult_part(psi_y_aux[i][j], b_in_y, psi_y_out, M/2, 3*M/4, G, M);

	for(int g = 0 ; g < G ; g++){
		for(int k = M/2 ; k < 3*M/4 ; k++){
			index = k + g*M;

			psi_x[i][j][index] = psi_x_out[index]   + (b_aux_x[i][j][index] - b1_x[i][j][index]);
			psi_y[i+1][j][index] = psi_y_out[index] + (b_aux_y[i][j][index] - b1_y[i][j][index]);
		}
	}
}

void sweep_3(const int &M, 
			 const int &G, 
			 int &nxf, 
			 int &nyf, 
			 vector<vector<vector<double>>> &l_x,
	 		 vector<vector<vector<double>>> &l_y, 
	 		 vector<double> &mi, 
	 		 vector<double> &n, 
	 		 vector<double> &w,
			 vector<vector<vector<double>>> &psi_x, 
			 vector<vector<vector<double>>> &psi_y, 
			 vector<vector<vector<double>>> &q,
			 vector<double> &b_in_x, 
			 vector<double> &b_in_y, 
			 vector<vector<double>> &h_x, 
			 vector<vector<double>> &h_y,
			 vector<double> &psi_x_out, 
			 vector<double> &psi_y_out,
			 vector<vector<vector<double>>> &c1_x, 
			 vector<vector<vector<double>>> &c2_x, 
			 vector<vector<vector<double>>> &c1_y, 
			 vector<vector<vector<double>>> &c2_y, 
			 vector<vector<vector<double>>> &b0_x, 
			 vector<vector<vector<double>>> &b1_x, 
			 vector<vector<vector<double>>> &b2_x, 
			 vector<vector<vector<double>>> &b0_y, 
			 vector<vector<vector<double>>> &b1_y, 
			 vector<vector<vector<double>>> &b2_y,
			 vector<vector<int>> &zcn, 
			 vector<vector<vector<double>>> &A_p_inv_z,
			 vector<vector<vector<double>>> &mhx, 
			 vector<vector<vector<double>>> &nhy, 
			 vector<vector<vector<double>>> &b_aux_x, 
			 vector<vector<vector<double>>> &b_aux_y,
			 vector<vector<vector<vector<double>>>> &psi_x_aux,
			 vector<vector<vector<vector<double>>>> &psi_y_aux){

	int index;

	quad_x_ne(nxf-1, 0, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_x);
	quad_y_ne(nxf-1, 0, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_y);

	get_psi_out_3(nxf-1, 0, M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y,
				q, b_in_x, b_in_y, h_x, h_y, psi_x_out, psi_y_out,
				b1_x, b1_y, zcn, mhx, nhy,
				b_aux_x, b_aux_y, psi_x_aux, psi_y_aux);

	for(int j = nxf-2 ; j >= 1 ; j--){
		quad_x(j, 0, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_x);
		quad_y_n(j, 0, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_y);

		get_psi_out_3(j, 0, M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y,
					q, b_in_x, b_in_y, h_x, h_y, psi_x_out, psi_y_out,
					b1_x, b1_y, zcn, mhx, nhy,
					b_aux_x, b_aux_y, psi_x_aux, psi_y_aux);
	}

	quad_x_nw(0, 0, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_x);
	quad_y_nw(0, 0, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_y);

	get_psi_out_3(0, 0, M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y,
				q, b_in_x, b_in_y, h_x, h_y, psi_x_out, psi_y_out,
				b1_x, b1_y, zcn, mhx, nhy,
				b_aux_x, b_aux_y, psi_x_aux, psi_y_aux);

	for(int i = 1 ; i < nyf-1 ; i++){
		
		quad_x_e(nxf-1, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_x);
		quad_y(nxf-1, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_y);

		get_psi_out_3(nxf-1, i, M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y, q, b_in_x, b_in_y, h_x, h_y, psi_x_out, psi_y_out,
				b1_x, b1_y, zcn, mhx, nhy, b_aux_x, b_aux_y, psi_x_aux, psi_y_aux);

		for(int j = nxf-2 ; j >= 1 ; j--){
			quad_x(j, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_x);
			quad_y(j, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_y);

			get_psi_out_3(j, i, M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y,	q, b_in_x, b_in_y, h_x, h_y, psi_x_out, psi_y_out,
					b1_x, b1_y, zcn, mhx, nhy, b_aux_x, b_aux_y, psi_x_aux, psi_y_aux);
		}

		quad_x_w(0, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_x);
		quad_y(0, i, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_y);

		get_psi_out_3(0, i, M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y,
				q, b_in_x, b_in_y, h_x, h_y, psi_x_out, psi_y_out,
				b1_x, b1_y, zcn, mhx, nhy,
				b_aux_x, b_aux_y, psi_x_aux, psi_y_aux);
	}

	quad_x_se(nxf-1, nyf-1, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_x);
	quad_y_se(nxf-1, nyf-1, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_y);

	get_psi_out_3(nxf-1, nyf-1, M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y,
				q, b_in_x, b_in_y, h_x, h_y, psi_x_out, psi_y_out,
				b1_x, b1_y, zcn, mhx, nhy,
				b_aux_x, b_aux_y, psi_x_aux, psi_y_aux);

	for(int j = nxf-2 ; j >= 1 ; j--){
		quad_x(j, nyf-1, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_x);
		quad_y_s(j, nyf-1, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_y);

		get_psi_out_3(j, nyf-1, M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y,
					q, b_in_x, b_in_y, h_x, h_y, psi_x_out, psi_y_out,
					b1_x, b1_y, zcn, mhx, nhy,
					b_aux_x, b_aux_y, psi_x_aux, psi_y_aux);
	}

	quad_x_sw(0, nyf-1, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_y, q, h_x, h_y, c1_x, c2_x, b0_x, b1_x, b2_x, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_x);
	quad_y_sw(0, nyf-1, M, G, nxf, nyf, mi, n, zcn, A_p_inv_z, l_x, q, h_x, h_y, c1_y, c2_y, b0_y, b1_y, b2_y, psi_x, psi_y, M/4, M/2, mhx, nhy, b_aux_y);

	get_psi_out_3(0, nyf-1, M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y,
				q, b_in_x, b_in_y, h_x, h_y, psi_x_out, psi_y_out,
				b1_x, b1_y, zcn, mhx, nhy,
				b_aux_x, b_aux_y, psi_x_aux, psi_y_aux);
	
}