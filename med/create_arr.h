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

	int R, Z, G, L, nxf, a = 0, iter = 1;

	double lb, rb, dev = 0, sum = 0;

	std::vector<int> nodes_R, zone_config, nodes;

	std::vector<double> mi, w, h_R, h, pos_arr, avg, alpha_temp, D, D_eig, DI_eig, q_p, psi_p_temp;

	std::vector<std::vector<double>> q_R, sigmat_R, sigmat_z, sigmat, psi, psi_m, phi, phi_old, q, SS, 
			alpha, A, ni, ni_i, V_temp, V_eig, VI_eig, temp, A_p_temp, psi_p;

	std::vector<std::vector<std::vector<double>>> sigmas0_R, sigmas1_R, sigmas0, sigmas1, sigmas0_z, sigmas1_z, V, VI, A_p, Am, Am_inv;