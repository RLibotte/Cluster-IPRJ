#include <iostream>
#include <vector>

tuple<vector<vector <double>>, vector<vector <double>>> getLeak(const int M, const int size_x, const int size_y, const int nxf, const int nyf, const int nx, const int ny,
								vector<vector<double>> h_x_nodes, vector<vector<double>> h_y_nodes, vector<vector<vector <double>>> psi_x, 
								vector<vector<vector <double>>> psi_y, vector<double> mi, vector<double> n, vector<double> w, vector<vector<int>> nodes_x,
								vector<vector<int>> nodes_y){

	vector<vector <double>> J_x (nyf, vector<double> (2,0)); vector<vector <double>> J_y (2, vector<double> (nxf,0));
	vector<vector <double>> J_x_out (size_y, vector<double> (2,0)); vector<vector <double>> J_y_out (2, vector<double> (size_x,0));

	for(int ky = 0 ; ky < size_y ; ky++){
		for(int i = 0 ; i < nodes_y[ky][0] ; i++){
			J_x[i][0] = 0;
			for(int j = 0 ; j < M/4 ; j++){
				J_x[i][0] = J_x[i][0] - mi[j + M/4] * psi_x[i][0][j + M/4] * h_y_nodes[i][0] * w[j + M/4];
				J_x[i][0] = J_x[i][0] - mi[j + M/2] * psi_x[i][0][j + M/2] * h_y_nodes[i][0] * w[j + M/2];
			}
			J_x_out[ky][0] = J_x_out[ky][0] + 0.5 * J_x[i][0];
		}
	}

	for(int ky = 0 ; ky < size_y ; ky++){
		for(int i = 0 ; i < nodes_y[ky][nodes_y.size() - 1] ; i++){
			J_x[i][1] = 0;
			for(int j = 0 ; j < M/4 ; j++){
				J_x[i][1] = J_x[i][1] + mi[j] * psi_x[i][nxf][j] * h_y_nodes[i][nx - 1] * w[j];
				J_x[i][1] = J_x[i][1] + mi[j + 3*M/4] * psi_x[i][nxf][j + 3*M/4] * h_y_nodes[i][nx - 1] * w[j + 3*M/4];
			}
			J_x_out[ky][1] = J_x_out[ky][1] + 0.5 * J_x[i][1];
		}	}


	for(int kx = 0 ; kx < size_x ; kx++){
		for(int i = 0 ; i < nodes_x[0][kx] ; i++){
			J_y[0][i] = 0;
			for(int j = 0 ; j < M/4 ; j++){
				J_y[0][i] = J_y[0][i] + n[j] * psi_y[0][i][j] * h_x_nodes[0][i] * w[j];
				J_y[0][i] = J_y[0][i] + n[j + M/4] * psi_y[0][i][j + M/4] * h_x_nodes[0][i] * w[j + M/4];
			}
			J_y_out[0][kx] = J_y_out[0][kx] + 0.5 * J_y[0][i];
		}
	}

	for(int kx = 0 ; kx < size_x ; kx++){
		for(int i = 0 ; i < nodes_x[nodes_y.size() - 1][kx] ; i++){
			J_y[1][i] = 0;
			for(int j = 0 ; j < M/4 ; j++){
				J_y[1][i] = J_y[1][i] - n[j + M/2] * psi_y[nyf][i][j + M/2] * h_x_nodes[ny - 1][i] * w[j + M/2];
				J_y[1][i] = J_y[1][i] - n[j + 3*M/4] * psi_y[nyf][i][j + 3*M/4] * h_x_nodes[ny - 1][i] * w[j + 3*M/4];
			}
			J_y_out[1][kx] = J_y_out[1][kx] + 0.5 * J_y[1][i];
		}
		cout << endl;
	}
	return {J_x_out, J_y_out};
}