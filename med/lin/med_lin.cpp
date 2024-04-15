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
#include <thread>

using namespace std;

#include "../../print/printarray.cpp"
#include "../../print/print_output.cpp"
#include "../../print/print_out.cpp"

#include "../../functions/get_model.cpp"

#include "../../eig/eig.cpp"
#include "../../functions/solve.cpp"

#include "../../leakage/get_lin_l.cpp"
#include "../../leakage/get_transverse_leakage.cpp"

#include "./get_matrix.cpp"
#include "./sweep_1.cpp"
#include "./sweep_2.cpp"
#include "./sweep_3.cpp"
#include "./sweep_4.cpp"
#include "./get_psi_in.cpp"

#include "../../functions/get_ref.cpp"
#include "../../functions/get_scalar_flux.cpp"
#include "../../functions/get_psi_p.cpp"
#include "../../functions/get_psi_b.cpp"
#include "../../functions/get_arrays.cpp"
#include "../../functions/get_eig.cpp"
#include "../../functions/getLQn.cpp"
#include "../../functions/getDev.cpp"
#include "../../functions/getS.cpp"
#include "../../functions/getLeakage_psi.cpp"
#include "../../functions/get_abs.cpp"

double med_lin(string model, const int N, const int M, const double max_dev, const int wx, const int wy){

	clock_t tStart = clock();

	auto [mi, n, w, order] = getLQn(N);

	auto [Z, x_size, y_size, ref_o, ref_s, ref_l, ref_n, nodes_x, nodes_y,  sigmat_z, sigmas0_z, sigmas1_z, 
	  q_n, h_x_n, h_y_n, G, zone_config] = get_model(model, wx, wy);

	auto [nxf, nyf, q, h_x, h_y, zcn, z, index, dev_x, MG] = 
		 get_arrays(M, G, x_size, y_size, q_n, nodes_x, nodes_y, h_x_n, h_y_n, zone_config);
		 
	auto [D_x, D_y, V_x, V_y, A_p_inv] =	
		 get_eig(M, G, Z, x_size, y_size, nxf, nyf, zone_config, mi, n, w, sigmat_z, sigmas0_z, sigmas1_z, nodes_x, nodes_y);
		 
	auto [b_in_x, b_in_y, psi_x_p, psi_y_p, l_x_temp, l_y_temp, psi_x_out, psi_y_out,
		  psi_x, psi_y, psi_b_x, psi_b_y, alpha_x, alpha_y, l_x, l_y, scalar_flux_x, 
		  scalar_flux_x_out, scalar_flux_y, scalar_flux_y_out, psi_x_p_out, psi_y_p_out, 
		  psi_x_old, psi_y_old,  
		  c1_x, b0_x, b1_x, c1_y, b0_y, b1_y, b_aux_x, b_aux_y, mhx, nhy, A_p_lu, 
		  psi_in_x, psi_in_y, psi_in_x_inv, psi_in_y_inv] = 
		  get_matrix(M, G, nxf, nyf, x_size, y_size, D_x, D_y, V_x, V_y, h_x, h_y, h_x_n, h_y_n, nodes_x, nodes_y, A_p_inv, zone_config, zcn, mi, n);

	get_psi_in_inv(nxf, nyf, M, G,  D_x, D_y, V_x, V_y, h_x, h_y, psi_in_x_inv, psi_in_y_inv, zone_config, zcn);
	get_psi_in(nxf, nyf, M, G,  D_x, D_y, V_x, V_y, h_x, h_y, psi_in_x, psi_in_y, zone_config, zcn);

	do{ ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		get_transverse_leakage(nxf, nyf, M, G, 3*M/4, M, psi_x, psi_y, mi, n, h_x, h_y, l_x, l_y);

		sweep_1(M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y, psi_x_old, psi_y_old, l_x_temp, l_y_temp, 
				q, b_in_x, b_in_y, alpha_x, alpha_y, psi_x_p, psi_y_p,
		 		psi_b_x, psi_b_y, h_x, h_y, psi_in_x, psi_in_y, psi_in_x_inv, psi_in_y_inv, 
		 		psi_x_p_out, psi_y_p_out, psi_x_out, psi_y_out,
				c1_x, c1_y, b0_x, b1_x, b0_y, b1_y, zcn, A_p_inv, mhx, nhy);

		if(ref_l == 1) get_ref_l(M, G, nxf, nyf, psi_x);
		if(ref_n == 1) get_ref_n(M, G, nxf, nyf, psi_y);

		get_transverse_leakage(nxf, nyf, M, G, 0, M/4, psi_x, psi_y, mi, n, h_x, h_y, l_x, l_y);

		sweep_2(M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y, psi_x_old, psi_y_old, l_x_temp, l_y_temp,  
		 		q, b_in_x, b_in_y, alpha_x, alpha_y, psi_x_p, psi_y_p,
		 		psi_b_x, psi_b_y, h_x, h_y, psi_in_x, psi_in_y, psi_in_x_inv, psi_in_y_inv, 
		 		psi_x_p_out, psi_y_p_out, psi_x_out, psi_y_out,
				c1_x, c1_y, b0_x, b1_x, b0_y, b1_y, zcn, A_p_inv, mhx, nhy);

		if(ref_o == 1) get_ref_o(M, G, nxf, nyf, psi_x);
		if(ref_n == 1) get_ref_n(M, G, nxf, nyf, psi_y);

		get_transverse_leakage(nxf, nyf, M, G, M/4, M/2, psi_x, psi_y, mi, n, h_x, h_y, l_x, l_y);

		sweep_3(M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y, psi_x_old, psi_y_old, l_x_temp, l_y_temp,
		 		q, b_in_x, b_in_y, alpha_x, alpha_y, psi_x_p, psi_y_p,
		 		psi_b_x, psi_b_y, h_x, h_y, psi_in_x, psi_in_y, psi_in_x_inv, psi_in_y_inv, 
		 		psi_x_p_out, psi_y_p_out, psi_x_out, psi_y_out,
				c1_x, c1_y, b0_x, b1_x, b0_y, b1_y, zcn, A_p_inv, mhx, nhy);

		if(ref_o == 1) get_ref_o(M, G, nxf, nyf, psi_x);
		if(ref_s == 1) get_ref_s(M, G, nxf, nyf, psi_y);

		get_transverse_leakage(nxf, nyf, M, G, M/2, 3*M/4, psi_x, psi_y, mi, n, h_x, h_y, l_x, l_y);

		sweep_4(M, G, nxf, nyf, l_x, l_y, mi, n, w, psi_x, psi_y, psi_x_old, psi_y_old, l_x_temp, l_y_temp,
		 		q, b_in_x, b_in_y, alpha_x, alpha_y, psi_x_p, psi_y_p,
		 		psi_b_x, psi_b_y, h_x, h_y, psi_in_x, psi_in_y, psi_in_x_inv, psi_in_y_inv, 
		 		psi_x_p_out, psi_y_p_out, psi_x_out, psi_y_out,
				c1_x, c1_y, b0_x, b1_x, b0_y, b1_y, zcn, A_p_inv, mhx, nhy);


		if(ref_l == 1) get_ref_l(M, G, nxf, nyf, psi_x);
		if(ref_s == 1) get_ref_s(M, G, nxf, nyf, psi_y);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// get_psi_p(nxf, nyf, G, M, l_x, l_y,	l_x_temp, l_y_temp, q, psi_x_p,	psi_y_p, A_p_lu, b_in_x, b_in_y, alpha_x, alpha_y, 
		// 		psi_in_x_inv, psi_in_y_inv, psi_x_p_out, psi_y_p_out, psi_x, psi_y);
		// get_psi_b(nxf, nyf, M, G, psi_b_x, psi_b_y, D_x, D_y, alpha_x, alpha_y, V_x, V_y, h_x, h_y, psi_x_p_out, psi_y_p_out, zcn);
		// get_scalar_flux(M, G, nxf, nyf, scalar_flux_x, scalar_flux_x_out, psi_b_x, nodes_x, nodes_y, w, x_size, y_size);
		// dev_x = getdev(x_size, y_size, G, scalar_flux_x_out, scalar_flux_x_old); //////////////////////////////////////////		
		// scalar_flux_x_old = scalar_flux_x_out;
		// print_flux(scalar_flux_x_out, "scalar_flux_x");

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		dev_x = getdev(nxf, nyf, MG, psi_x, psi_x_old); //////////////////////////////////////////		
		psi_x_old = psi_x;

		// z++, printf("%d --> %.6e \n", z, dev_x);

		// for(int i = 0 ; i < nyf ; i++){
		// 	for(int j = 0 ; j < nxf + 1 ; j++){
		// 		for(int g = 0 ; g < G ; g++){
		// 		psi_x_avg[i][j][g] = 0;
		// 			for(int l = 0 ; l < M ; l++){
		// 				psi_x_avg[i][j][g] += psi_x[i][j][l + M*g] * w[l];
		// 			}
		// 			// printf("%.6e ", psi_x_avg[i][j][g]);
		// 		}
		// 		// printf("\n");
		// 	}
		// }

		// print_flux(scalar_flux_x_out, "scalar_flux_y");

	} while(dev_x > max_dev);
	// } while(z < 100);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	get_psi_p(nxf, nyf, G, M, l_x, l_y,	l_x_temp, l_y_temp, q, psi_x_p,	psi_y_p, A_p_lu, b_in_x, b_in_y, alpha_x, alpha_y, 
		psi_in_x_inv, psi_in_y_inv, psi_x_p_out, psi_y_p_out, psi_x, psi_y);
	
	get_psi_b(nxf, nyf, M, G, psi_b_x, psi_b_y, D_x, D_y, alpha_x, alpha_y, V_x, V_y, h_x, h_y, psi_x_p_out, psi_y_p_out, zcn);
	
	get_scalar_flux(M, G, nxf, nyf, scalar_flux_x, scalar_flux_x_out, psi_b_x, nodes_x, nodes_y, w, x_size, y_size);
	get_scalar_flux(M, G, nxf, nyf, scalar_flux_y, scalar_flux_y_out, psi_b_y, nodes_x, nodes_y, w, x_size, y_size);

	for(int i = 0 ; i < y_size ; i++){
		for(int j = 0 ; j < x_size ; j++){
			for(int g = 0 ; g < G ; ++g){
				scalar_flux_x_out[i][j][g] = (scalar_flux_x_out[i][j][g] + scalar_flux_y_out[i][j][g])/2;
			}
		}
	}

	// print_flux(scalar_flux_x_out, "scalar_flux"); 

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// print_flux(psi_x_avg, "psi_x_avg");

	// print_flux(psi_x, "psi_x");
	// print_flux(psi_y, "psi_y");

	// print_flux(scalar_flux_x_out, "scalar_flux_x");

    auto [t_abs_x_out, t_abs_x] = 
    	get_abs(G, x_size, y_size, nxf, nyf, sigmat_z, sigmas0_z, h_x_n, h_y_n, h_x, h_y, scalar_flux_x, scalar_flux_x_out, zone_config, zcn);

	// print_flux(t_abs_x_out, "t_abs_x_out");

	// auto [J_x, J_y] = getLeak_psi(M, x_size, y_size, nxf, nyf, G, h_x, h_y, psi_x, psi_y, mi, n, w, nodes_x, nodes_y);
		// print_leakage_out(M, G, x_size, y_size, J_x, J_y);

	//////////////////////////////////////////////////////////////////////

	// // model 01

	// double ref_s4 = 1.303422e+01; //7.11992334e-01; ;
    // double ref_s4 = 1.28244563e+01; // S2

	
	// double ref_s4 = 7.30136309e-01 ;// S8 ------------

	// double jd = 3.482862e+01;

	// // double ref_s8 = 1.30702395e+01;
	// // double jd = 3.46474300e+01;



	// printf("%d $\\times$ %d & %.8e (%.6f $\\%$) & %.8e (%.6f $\\%$) \\\\ \n", wx, wx, scalar_flux_x_out[0][0][0], fabs(ref_s4 - scalar_flux_x_out[0][0][0])/ref_s4 * 100, J_x[0][1], fabs(jd - J_x[0][1])/jd * 100);

	// printf("%d x %d --- %.10e --- dev =  %.8f \n", wx, wy, scalar_flux_x_out[0][0][0], fabs(ref_s4 - scalar_flux_x_out[0][0][0])/ref_s4 * 100);


	// printf("%.8f,\n", fabs(ref_s4 - scalar_flux_x_out[0][0][0])/ref_s4 * 100);

	//////////////////////////////////////////////////////////////////////

	// // model 20

	// double ref_s6 = 4.347705e-02;
	// double ref_s16 = 4.355931e-02;


	// double ref_s16 = 4.355931e-02;

	// printf("wx =  %d --- %.6e --- deviation = %.6f % \n", wx, t_abs_x_out[1][2][0], fabs(ref_s16 - t_abs_x_out[1][2][0])/ref_s16 * 100);

	///////////// TESE DANY

	// double ref_s6_d1 = 1.716E+00;
	// double ref_s6_d2 = 1.246E-02;

	// printf("%d $\\times$ %d & %.8e (%.6f $\\%$) & %.8e (%.6f $\\%$) \\\\ \n", wx*7, wx*8, scalar_flux_x_out[3][2][0], fabs(ref_s6_d1 - scalar_flux_x_out[3][2][0])/ref_s6_d1 * 100, scalar_flux_x_out[1][2][0], fabs(ref_s6_d2 - scalar_flux_x_out[1][2][0])/ref_s6_d2 * 100);
	// printf("wx =  %d --- %.6e --- deviation = %.6f % \n", wx, t_abs_x_out[1][2][0], fabs(ref_s16 - t_abs_x_out[1][2][0])/ref_s16 * 100);


	//////////////////////////////////////////////////////////////////////

	// model 26 

	

	// vector<double> ref (G);

	// S4

	// ref[0] = 4.116687e-04;
	// ref[1] = 6.543834e-04;

 	// S8

	// ref[0] = 4.386200e-04;
	// ref[1] = 7.052616e-04;
 	

	// printf("$%d \\times %d$ ", nxf, nyf);

	// double temp_f = 0;
	// for(int g = 0 ; g < G ; g++){
	// 	temp_f = 0;
	// 	for(int i = 0 ; i < y_size ; i++){
	// 		temp_f = temp_f +  J_x[i][g + G];
	// 	}

	// 	printf(" & %.6e (%.6f \\%) ", temp_f, fabs(ref[g]-temp_f)/ref[g] * 100);
	// }
	// printf(" \\\\ \n");

	//////////////////////////////////////////////////////////////////////

	// printf("%d x %d & %.6e & %.6e & %.6e \\\\ \n", nxf, nyf, scalar_flux_x_out[0][3][0], scalar_flux_x_out[0][3][4], scalar_flux_x_out[0][3][9]);

	// print_output(N, x_size, y_size, nxf, nyf, 1, 1, G, scalar_flux_x, t_abs_x, scalar_flux_x, t_abs_x_out, "med", model);

	// printf("%.6f \n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

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

	// int wx;
	// cout << "wx = ";
	// cin >> wx;
	
	// const int wx = 2;
	// const int wy = wx;

	string model = "01";

	// med_2d(model, N, M, max_dev, wx, wy);

	vector<int> wx_vec = {50};

	double time;
	int iter = 5;
	for(int i = 0 ; i < wx_vec.size() ; i++){
		printf("wx = %d \n\n", wx_vec[i]);
		time = 0;
		for(int j = 0 ; j < iter ; j++){
			time += med_lin(model, N, M, max_dev, wx_vec[i], wx_vec[i]);
		}
		printf(" %6f \n\n", time/iter);
	}

	// for(int i = 0 ; i < wx_vec.size() ; i++){
	// 	med_2d(model, N, M, max_dev, wx_vec[i], wx_vec[i]);
	// }

	printf("Time taken: %.6f s\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	return 0;
}
