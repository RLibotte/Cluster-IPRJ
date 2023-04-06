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

	int R, Z, G, L, nxf, a = 0;

	double lb, rb, dev = 0;

	std::vector<int> nodes_R, zone_config, nodes;

	std::vector<double> mi, w, h_R, h, pos_arr;

	std::vector<std::vector<double>> q_R, sigmat_R, sigmat, psi, psi_m, phi, phi_old, sigmat_z, q;

	std::vector<std::vector<std::vector<double>>> sigmas0_R, sigmas1_R, sigmas0, sigmas1, SS, sigmas0_z, sigmas1_z;