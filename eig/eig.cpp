#include <iostream>
#include <string>
#include <vector>

#include "./Alglib/linalg.h"
#include "./Alglib/solvers.h"
#include "./Alglib/ap.h"

using namespace std;

tuple<vector<vector<double>>, vector<double>> 
    eigen(vector<vector<double>> &arr, const int NG) {

    vector<vector <double>> V(NG , vector<double> (NG)); 
    vector<double> D(NG);
    
    int VNeeded = 3;
    bool Conv;
    double *VV; 
    alglib::real_2d_array AA, VR, VL;
    alglib::real_1d_array WR, WI;

    VV = new double[NG * NG];
    VR.setlength(NG, NG);
    VL.setlength(NG, NG);
    WR.setlength(NG);
    WI.setlength(NG);

    int k = 0;
    for (int i = 0; i < NG; i++){
        for (int j = 0; j < NG; j++) {
            VV[j + NG * i] = arr[i][j];
        }
    }

    AA.setcontent(NG, NG, VV);
    Conv = alglib::rmatrixevd(AA, NG, VNeeded, WR, WI, VL, VR);

    for(int i = 0 ; i < NG ; i++){
        for(int j = 0 ; j < NG ; j++){
            V[i][j] = VR[i][j];
        }
        D[i] = WR[i];
    }

    free(VV);
    
    return {V, D};
}



tuple<vector<vector<double>>, vector<double>, vector<double>> 
    eigen_c(vector<vector<double>> &arr, const int NG) {

    vector<vector <double>> V(NG , vector<double> (NG)); 
    vector<double> D(NG);
    vector<double> D_C(NG);
    
    int VNeeded = 3;
    bool Conv;
    double *VV; 
    alglib::real_2d_array AA, VR, VL;
    alglib::real_1d_array WR, WI;

    VV = new double[NG * NG];

    VR.setlength(NG, NG);
    VL.setlength(NG, NG);
    WR.setlength(NG);
    WI.setlength(NG);

    for (int i = 0; i < NG; i++){
        for (int j = 0; j < NG; j++) {
            VV[j + i*NG] = arr[i][j];
        }
    }

    AA.setcontent(NG, NG, VV);
    Conv = alglib::rmatrixevd(AA, NG, VNeeded, WR, WI, VL, VR);

    for(int i = 0 ; i < NG ; i++){
        for(int j = 0 ; j < NG ; j++){
            V[i][j] = VR[i][j];

            // printf("%.6f ", VL[i][j]);
        }
        D[i] = WR[i];
        D_C[i] = WI[i];
        // princtf("\n");
        if(WI[i] != 0) printf("WI -> %.6e \n", WI[i]);
    }
    // printf("\n");

    free(VV);
    
    return {V, D, D_C};
}
