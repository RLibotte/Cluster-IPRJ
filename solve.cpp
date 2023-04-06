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

vector<vector <double>> inverse(vector<vector <double>> arr, const int NG) {

    vector<vector <double>> ans(NG , vector<double> (NG)); 
    alglib::real_2d_array A;
    alglib::ae_int_t info;
    alglib::matinvreport rep;

    A.setlength(NG, NG);

    double* MA = new double[NG*NG];

    int k = 0;
    for (int i = 0; i < NG; i++){
        for (int j = 0; j < NG; j++) {
            MA[k] = arr[i][j];
            k++;
        }
    }
    
    A.setcontent(NG, NG, MA);

    free(MA);

    alglib::rmatrixinverse(A, info, rep);

    for (int i = 0; i < NG; i++) {
        for(int j = 0 ; j < NG ; j++){
            ans[i][j] = A[i][j];
        }
    }

    // print2dArray(ans);

    return ans;
}

vector<vector <double>> luinverse(vector<vector <double>> arr, const int NG) {

    vector<vector <double>> ans(NG , vector<double> (NG));
    alglib::real_2d_array lua;
    alglib::ae_int_t info, n;
    alglib::matinvreport rep;
    alglib::integer_1d_array p;

    lua.setlength(NG, NG);

    double* MA = new double[NG * NG];

    int k = 0;
    for (int i = 0; i < NG; i++){
        for (int j = 0; j < NG; j++) {
            MA[k] = arr[i][j];
            k++;
        }
    }
    
    lua.setcontent(NG, NG, MA);

    free(MA);

    alglib::rmatrixlu(lua, NG, NG, p);
    alglib::rmatrixluinverse(lua, p, info, rep);

    for (int i = 0; i < NG; i++) {
        for(int j = 0 ; j < NG ; j++){
            ans[i][j] = lua[i][j];
        }
    }

    return ans;
}

vector<double> mult(vector<vector <double>> mat, vector<double> arr){

    const int size = mat.size();

    vector<double> ans(size,0);

    double sum;

    for(int i = 0 ; i < size ; i++){
        sum = 0;
    //    #pragma omp parallel for reduction(+:sum)
        for(int j = 0 ; j < size ; j++){
            sum += mat[i][j] * arr[j];
        }
        ans[i] = sum;
    }

    return ans;

}

vector<vector<double>> mult_mat(vector<vector <double>> mat_1, vector<vector<double>> mat_2){

    const int size = mat_1.size();
    vector<vector<double>> ans (size, vector<double> (size));

    for(int i = 0 ; i < size ; i++){
        for(int j = 0 ; j < size ; j++){
            for(int k = 0 ; k < size ; k++){
                ans[i][j] += mat_1[i][k] * mat_2[k][j];
            }
        }
    }

    return ans;

}