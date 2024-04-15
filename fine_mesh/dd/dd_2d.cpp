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

#include "./get_matrix.cpp"
#include "./get_psi_b.cpp"
#include "./get_s.cpp"

#include "./sweep_1.cpp"
#include "./sweep_2.cpp"
#include "./sweep_3.cpp"
#include "./sweep_4.cpp"

#include "../../functions/get_ref.cpp"
#include "../../functions/get_scalar_flux.cpp"
#include "../../functions/get_arrays.cpp"
#include "../../functions/getLQn.cpp"
#include "../../functions/getDev.cpp"
#include "../../functions/getS.cpp"
#include "../../functions/getLeakage_psi.cpp"
#include "../../functions/get_abs.cpp"

void dd_2d(string model, const int N, const int M, const double max_dev, const int wx, const int wy){

	auto [mi, n, w, order] = getLQn(N);

	auto [Z, x_size, y_size, ref_o, ref_s, ref_l, ref_n, nodes_x, nodes_y,  sigmat_z, sigmas0_z, sigmas1_z, 
	  q_n, h_x_n, h_y_n, G, zone_config] = get_model(model, wx, wy);

	auto [nxf, nyf, q, h_x, h_y, zcn, z, index, dev_x, MG] = 
		 get_arrays(M, G, x_size, y_size, q_n, nodes_x, nodes_y, h_x_n, h_y_n, zone_config);

	auto [psi_x_b, psi_y_b, scalar_flux_x, scalar_flux_x_old, scalar_flux_x_out, SS_x, SS_y, SS_x_1, SS_y_1, psi_x, psi_x_old, psi_y] =
		 get_matrix(x_size, y_size, nxf, nyf, M, G);

 	do{

	 	get_psi_b(M*G, nxf, nyf, psi_x, psi_y, psi_x_b, psi_y_b);
		get_s(M, G, nxf, nyf, w, mi, n, sigmas0_z, sigmas1_z, psi_x_b, psi_y_b, SS_x, SS_y, SS_x_1, SS_y_1, zcn);

		sweep_1(M, G, nxf, nyf, mi, n, w, psi_x_b, psi_y_b, psi_x, psi_y, h_x, h_y, SS_x, SS_y, SS_x_1, SS_y_1, q, sigmat_z, sigmas1_z, zcn);
		sweep_2(M, G, nxf, nyf, mi, n, w, psi_x_b, psi_y_b, psi_x, psi_y, h_x, h_y, SS_x, SS_y, SS_x_1, SS_y_1, q, sigmat_z, sigmas1_z, zcn);
		sweep_3(M, G, nxf, nyf, mi, n, w, psi_x_b, psi_y_b, psi_x, psi_y, h_x, h_y, SS_x, SS_y, SS_x_1, SS_y_1, q, sigmat_z, sigmas1_z, zcn);
		sweep_4(M, G, nxf, nyf, mi, n, w, psi_x_b, psi_y_b,	psi_x, psi_y, h_x, h_y, SS_x, SS_y, SS_x_1, SS_y_1, q, sigmat_z, sigmas1_z, zcn);

		if(ref_l == 1) get_ref_l(M, G, nxf, nyf, psi_x);
		if(ref_n == 1) get_ref_n(M, G, nxf, nyf, psi_y);
		if(ref_o == 1) get_ref_o(M, G, nxf, nyf, psi_x);
		if(ref_s == 1) get_ref_s(M, G, nxf, nyf, psi_y);

		z++;

		dev_x = getdev(nxf, nyf, MG, psi_x, psi_x_old);
		psi_x_old = psi_x;

		get_scalar_flux(M, G, nxf, nyf, scalar_flux_x, scalar_flux_x_out, psi_x_b, nodes_x, nodes_y, w, x_size, y_size);
		print_flux(scalar_flux_x_out, "scalar_flux_x_out");

 		printf("%d --> %.6e \n", z, dev_x);

 	} while(dev_x > max_dev);
 	// } while(z < 100);

	get_scalar_flux(M, G, nxf, nyf, scalar_flux_x, scalar_flux_x_out, psi_x_b, nodes_x, nodes_y, w, x_size, y_size);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// print_flux(scalar_flux_x_out, "scalar_flux_x_out");

 	// print_flux(psi_x, "psi_x");
 	// print_flux(psi_y, "psi_y");

	// print_flux(scalar_flux_x, "scalar_flux_x");
	// print_psi_out(M, G, wy, wx, nxf, nyf, psi_x, psi_y);
	
	////////////////////////////////////////////////////////////////////////////////////////

	// auto [J_x, J_y] = getLeak_psi(M, x_size, y_size, nxf, nyf, G, h_x, h_y, psi_x, psi_y, mi, n, w, nodes_x, nodes_y);
	// 	print_leakage_out(M, G, x_size, y_size, J_x, J_y);

	// auto [J_x, J_y] = getLeak_psi(M, x_size, y_size, nxf, nyf, G, h_x, h_y, psi_x, psi_y, mi, n, w, nodes_x, nodes_y);
	// 	print_leakage_out(M, G, x_size, y_size, J_x, J_y);

	// print2d(J_x, "J_x");

	////////////////////////////////////////////////////////////////////////////////////////

    // auto [t_abs_x_out, t_abs_x] = 
    // 	get_abs(G, x_size, y_size, nxf, nyf, sigmat_z, sigmas0_z, h_x_n, h_y_n, h_x, h_y, scalar_flux_x, scalar_flux_x_out, zone_config, zcn);

    ////////////////////////////////////////////////////////////////////////////////////////

	// print_flux(scalar_flux_x_out, "scalar_flux_x_out");
	// print_flux(t_abs_x_out, "t_abs_x_out");
	// print_flux(t_abs_x, "t_abs_x");

	// printf("%d x %d & %.6e & %.6e & %.6e \\\\ \n", nxf, nyf, scalar_flux_x_out[0][3][0], scalar_flux_x_out[0][3][4], scalar_flux_x_out[0][3][9]);

	// print_output(N, x_size, y_size, nxf, nyf, 1, 1, G, scalar_flux_x, t_abs_x, scalar_flux_x, t_abs_x_out, "dd", model);

    // system("python3 ../print_chart/print_chart.py");

}

int main(){

	clock_t tStart = clock();

    const int num_threads = 1;
	omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
    	printf("%d THREADS\n", omp_get_num_threads());
    }

	const int N = 8;
	const int M = N*(N+2)/2;

	double max_dev = pow(10, -7);

	const int wx = 100;
	const int wy = wx;

	string model = "19";

	dd_2d(model, N, M, max_dev, wx, wy);

	printf("Time taken: %.6f s\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	return 0;
}