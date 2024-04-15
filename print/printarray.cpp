#include <iostream>
#include <vector>

using namespace std;

void print1d(vector <double> arr, string name){
	cout << name << endl << endl;
	for(int i = 0 ; i < arr.size() ; i++){
		printf("%.6f \n", arr[i]);
	}
	printf("\n\n");
}

void print2d(vector<vector <double>> arr, string name){
	cout << name << endl << endl;
	for(int i = 0 ; i < arr.size() ; i++){
		for(int j = 0 ; j < arr[0].size() ; j++){
			printf("%.6f \t", arr[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
}

void print_flux(vector<vector<vector <double>>> arr, string name){
	cout << name << endl << endl;

	for(int i = 0 ; i < arr[0][0].size() ; i++){
		for(int j = 0 ; j < arr.size() ; j++){
			for(int k = 0 ; k < arr[0].size() ; k++){
				printf("%.8e \t", arr[j][k][i]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("\n\n");
}

void print3d(vector<vector<vector <double>>> arr, string name){
	cout << name << endl << endl;

	for(int i = 0 ; i < arr.size() ; i++){
		for(int j = 0 ; j < arr[0].size() ; j++){
			for(int k = 0 ; k < arr[0][0].size() ; k++){
				printf("%.6f \t", arr[i][j][k]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("\n\n");
}

void print4d(vector<vector<vector<vector<double>>>> arr, string name){
	cout << name << endl << endl;
	
	for(int i = 0 ; i < arr.size() ; i++){
		for(int j = 0 ; j < arr[0].size() ; j++){
			for(int k = 0 ; k < arr[0][0].size() ; k++){	
				for(int l = 0 ; l < arr[0][0][0].size() ; l++){
					printf("%.6f \t", arr[i][j][k][l]);
				}
				printf("\n");
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("\n\n");
}