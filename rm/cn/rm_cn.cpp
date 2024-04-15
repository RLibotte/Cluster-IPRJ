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

#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <tuple>
#include <chrono>
#include <omp.h>

using namespace std;

#include "../../print/printarray.cpp"
#include "../../print/print_output.cpp"
#include "../../print/print_out.cpp"

#include "../../functions/get_model.cpp"

#include "../../eig/eig.cpp"
#include "../../functions/solve.cpp"

#include "../../leakage/get_transverse_leakage.cpp"

#include "./get_matrix.cpp"

#include "./get_b_in_cn.cpp"
#include "./sweep_1.cpp"
#include "./sweep_2.cpp"
#include "./sweep_3.cpp"
#include "./sweep_4.cpp"
#include "./get_psi_b.cpp"

#include "../get_psi_in.cpp"

#include "../../functions/get_ref.cpp"
#include "../../functions/get_scalar_flux.cpp"
#include "../../functions/get_psi_p.cpp"
#include "../../functions/get_arrays.cpp"
#include "../../functions/get_eig.cpp"
#include "../../functions/getLQn.cpp"
#include "../../functions/getDev.cpp"
#include "../../functions/getS.cpp"
#include "../../functions/getLeakage_psi.cpp"
#include "../../functions/get_abs.cpp"

double rm_cn(string model, const int N, const int M, const double max_dev, const int wx, const int wy){

	clock_t tStart = clock();

	auto [mi, n, w, order] = getLQn(N);

	auto [Z, x_size, y_size, ref_o, ref_s, ref_l, ref_n, nodes_x, nodes_y,  sigmat_z, sigmas0_z, sigmas1_z, 
	  q_n, h_x_n, h_y_n, G, zone_config] = get_model(model, wx, wy);

	auto [nxf, nyf, q, h_x, h_y, zcn, z, index, dev_x, MG] = 
		 get_arrays(M, G, x_size, y_size, q_n, nodes_x, nodes_y, h_x_n, h_y_n, zone_config);

	auto [D_x, D_y, V_x, V_y, A_p_inv] =	
		 get_eig(M, G, Z, x_size, y_size, nxf, nyf, zone_config, mi, n, w, sigmat_z, sigmas0_z, sigmas1_z, nodes_x, nodes_y);

	auto [b_in_x, b_in_y, psi_x_p, psi_y_p, l_x_temp, l_y_temp, psi_x_out, psi_y_out,
			psi_x, psi_x_old, psi_y, psi_b_x, psi_b_y, alpha_x, alpha_y, l_x, l_y, l_x_old, scalar_flux_x, scalar_flux_x_out, psi_x_p_out,
			psi_y_p_out, mi_h, n_h,
			A_p_lu, psi_in_x, psi_in_y, psi_in_x_inv, psi_in_y_inv, psi_iter_x, psi_iter_y, psi_x_aux, psi_y_aux] = 
		  get_matrix(M, G, nxf, nyf, x_size, y_size, D_x, D_y, V_x, V_y, h_x, h_y, h_x_n, h_y_n, nodes_x, nodes_y, A_p_inv, zone_config, zcn, mi, n);

	get_psi_in_inv(nxf, nyf, M, G,  D_x, D_y, V_x, V_y, h_x, h_y, psi_in_x_inv, psi_in_y_inv, zone_config, zcn);
	get_psi_in(nxf, nyf, M, G,  D_x, D_y, V_x, V_y, h_x, h_y, psi_in_x, psi_in_y, zone_config, zcn);
	
	for(int i = 0 ; i < nyf ; i++){
		for(int j = 0 ; j < nxf ; j++){
			psi_iter_x[i][j] = mult_mat(psi_in_x[i][j], psi_in_x_inv[i][j]);
			psi_iter_y[i][j] = mult_mat(psi_in_y[i][j], psi_in_y_inv[i][j]);

			psi_x_aux[i][j] = add_eye(psi_iter_x[i][j]);
			psi_y_aux[i][j] = add_eye(psi_iter_y[i][j]);

			psi_x_aux[i][j] = mult_mat(psi_x_aux[i][j], A_p_lu[i][j]);
			psi_y_aux[i][j] = mult_mat(psi_y_aux[i][j], A_p_lu[i][j]);
		}
	}

	do{ ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		get_transverse_leakage(nxf, nyf, M, G, 3*M/4, M, psi_x, psi_y, mi, n, h_x, h_y, l_x, l_y);

		sweep_1(M, G, nxf, nyf, l_x, l_y, psi_x, psi_y, mi_h, n_h, l_x_temp, l_y_temp, 
				q, b_in_x, b_in_y, h_x, h_y, psi_iter_x, psi_iter_y, psi_x_aux, psi_y_aux);

		if(ref_l == 1) get_ref_l(M, G, nxf, nyf, psi_x);
		if(ref_n == 1) get_ref_n(M, G, nxf, nyf, psi_y);

		get_transverse_leakage(nxf, nyf, M, G, 0, M/4, psi_x, psi_y, mi, n, h_x, h_y, l_x, l_y);

		sweep_2(M, G, nxf, nyf, l_x, l_y, psi_x, psi_y, mi_h, n_h, l_x_temp, l_y_temp, 
				q, b_in_x, b_in_y, h_x, h_y, psi_iter_x, psi_iter_y, psi_x_aux, psi_y_aux);

		if(ref_o == 1) get_ref_o(M, G, nxf, nyf, psi_x);
		if(ref_n == 1) get_ref_n(M, G, nxf, nyf, psi_y);

		get_transverse_leakage(nxf, nyf, M, G, M/4, M/2, psi_x, psi_y, mi, n, h_x, h_y, l_x, l_y);

		sweep_3(M, G, nxf, nyf, l_x, l_y, psi_x, psi_y, mi_h, n_h, l_x_temp, l_y_temp, 
				q, b_in_x, b_in_y, h_x, h_y, psi_iter_x, psi_iter_y, psi_x_aux, psi_y_aux);

		if(ref_o == 1) get_ref_o(M, G, nxf, nyf, psi_x);
		if(ref_s == 1) get_ref_s(M, G, nxf, nyf, psi_y);

		get_transverse_leakage(nxf, nyf, M, G, M/2, 3*M/4, psi_x, psi_y, mi, n, h_x, h_y, l_x, l_y);

		sweep_4(M, G, nxf, nyf, l_x, l_y, psi_x, psi_y, mi_h, n_h, l_x_temp, l_y_temp, 
				q, b_in_x, b_in_y, h_x, h_y, psi_iter_x, psi_iter_y, psi_x_aux, psi_y_aux);

		if(ref_l == 1) get_ref_l(M, G, nxf, nyf, psi_x);
		if(ref_s == 1) get_ref_s(M, G, nxf, nyf, psi_y);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		dev_x = getdev(nxf, nyf, MG, psi_x, psi_x_old);
		psi_x_old = psi_x;

		// z++, printf("%d --> %.6e \n", z, dev_x);

	} while(dev_x > max_dev);
	// } while(z < 1);

	get_psi_p(nxf, nyf, G, M, l_x, l_y,	l_x_temp, l_y_temp, q, psi_x_p,	psi_y_p, A_p_lu, b_in_x, b_in_y, alpha_x, alpha_y, 
		psi_in_x_inv, psi_in_y_inv, psi_x_p_out, psi_y_p_out, psi_x, psi_y);

	get_psi_b(nxf, nyf, M, G, psi_b_x, psi_b_y, D_x, D_y, alpha_x, alpha_y, V_x, V_y, h_x, h_y, psi_x_p_out, psi_y_p_out, zcn);
	get_scalar_flux(M, G, nxf, nyf, scalar_flux_x, scalar_flux_x_out, psi_b_x, nodes_x, nodes_y, w, x_size, y_size);

	// print_flux(psi_x, "psi_x");
	// print_flux(scalar_flux_x_out, "scalar_flux");

    // auto [t_abs_x_out, t_abs_x] = 
    // 	get_abs(G, x_size, y_size, nxf, nyf, sigmat_z, sigmas0_z, h_x_n, h_y_n, h_x, h_y, scalar_flux_x, scalar_flux_x_out, zone_config, zcn);

	// print_flux(t_abs_x_out, "t_abs_x_out");

	auto [J_x, J_y] = getLeak_psi(M, x_size, y_size, nxf, nyf, G, h_x, h_y, psi_x, psi_y, mi, n, w, nodes_x, nodes_y);
		// print_leakage_out(M, G, x_size, y_size, J_x, J_y);

	// print_output(x_size, y_size, nxf, nyf, 17, 17, G, scalar_flux_x, t_abs_x, scalar_flux_x, t_abs_x_out, "med");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// model 1

	// double ref_s4 = 1.303405e+01; // S4
	// double jd = 3.482838e+01; // S4

	// double ref_s4 = 1.307023e+01; // S8
	// double jd = 3.464743e+01; // S8

	double ref_s4 = 1.308099e+01; // S16
	double jd = 3.459394e+01; // S16

	printf(" & %d $\\times$ %d & %.6e (%.6f $\\%$) & %.6e (%.6f $\\%$) \\\\", wx, wx, scalar_flux_x_out[0][0][0], fabs(ref_s4 - scalar_flux_x_out[0][0][0])/ref_s4 * 100, J_x[0][1], fabs(jd - J_x[0][1])/jd * 100);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// model 19

	// double ref_s6 = 2.72428480e-08; // S4
	// double ref_s6 = 3.88206368e-08; // S8
	// double ref_s6 = 4.355931e-02; // S16

	// printf(" & %d $\\times$ %d & %.6e (%.6f $\\%$) \\\\ \n", 10*wx, 10*wx, scalar_flux_x_out[0][4][0], fabs(ref_s6 - scalar_flux_x_out[0][4][0])/ref_s6 * 100);

	// printf("wx =  %d --- %.6e --- deviation = %.6f %% \n", wx, t_abs_x_out[1][2][0], fabs(ref_s6 - t_abs_x_out[1][2][0])/ref_s6 * 100);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// double ref_s6 = 1.24876510e-02; // S16
	// printf(" & %d $\\times$ %d & %.6e (%.6f $\\%$) \\\\ \n", 7*wx, 8*wx, scalar_flux_x_out[1][2][0], fabs(ref_s6 - scalar_flux_x_out[1][2][0])/ref_s6 * 100);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// model 26

	// vector<double> ref (G);

	// S4

	// ref[0] = 4.116687e-04;
	// ref[1] = 6.543834e-04;

 	// S8

	// ref[0] = 4.386200e-04;
	// ref[1] = 7.052616e-04;
 	

	// printf(" & $%d \\times %d$ ", nxf, nyf);

	// double temp_f = 0;
	// for(int g = 0 ; g < G ; g++){
	// 	temp_f = 0;
	// 	for(int i = 0 ; i < y_size ; i++){
	// 		temp_f = temp_f +  J_x[i][g + G];
	// 	}
	// 	printf(" & %.6e (%.6f \\%) ", temp_f, fabs(ref[g]-temp_f)/ref[g] * 100);
	// }
	// printf(" \\\\ \n");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// model 20

	// double ref_s6 = 1.24173742e-02; // S4
	// double ref_s6 = 1.24876510e-02; // S16

	// printf(" & %d $\\times$ %d & %.6e (%.6f $\\%$) \\\\ \n", 7*wx, 8*wx, scalar_flux_x_out[1][2][0], fabs(ref_s6 - scalar_flux_x_out[1][2][0])/ref_s6 * 100);

	// printf("%.6e", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	return (double)(clock() - tStart)/CLOCKS_PER_SEC;

}

int main(){

	clock_t tStart = clock();

    const int num_threads = 1;
	omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
    	printf("%d THREADS\n", omp_get_num_threads());
    }

	const int N = 16;
	const int M = N*(N+2)/2;

	double max_dev = pow(10, -7);

	string model = "26";

	const int wx = 50;
	const int wy = wx;

	rm_cn(model, N, M, max_dev, wx, wy);

	// vector<int> wx_vec = {1,2,3,4,5};

	// double time;
	// int iter = 10;
	// for(int i = 0 ; i < wx_vec.size() ; i++){
	// 	printf("wx = %d \n\n", wx_vec[i]);
	// 	time = 0;
	// 	for(int j = 0 ; j < iter ; j++){
	// 		time += rm_cn(model, N, M, max_dev, wx_vec[i], wx_vec[i]);
	// 	}
	// 	printf(" %6f \n\n", time/10);
	// }

	// for(int i = 0 ; i < wx_vec.size() ; i++){
	// 	rm_cn(model, N, M, max_dev, wx_vec[i], wx_vec[i]);
	// }

	printf("Time taken: %.6f s\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	return 0;
}
