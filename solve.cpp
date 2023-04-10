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

std::vector<std::vector <double>> inverse(std::vector<std::vector <double>> arr, const int NG) {

    std::vector<std::vector <double>> ans(NG , std::vector<double> (NG)); 
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

std::vector<std::vector <double>> luinverse(std::vector<std::vector <double>> arr, const int NG) {

    std::vector<std::vector <double>> ans(NG , std::vector<double> (NG));
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

std::vector<double> mult(std::vector<std::vector <double>> mat, std::vector<double> arr){

    const int size = mat.size();

    std::vector<double> ans(size,0);

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

std::vector<std::vector<double>> mult_mat(std::vector<std::vector <double>> mat_1, std::vector<std::vector<double>> mat_2){

    const int size = mat_1.size();
    std::vector<std::vector<double>> ans (size, std::vector<double> (size));

    for(int i = 0 ; i < size ; i++){
        for(int j = 0 ; j < size ; j++){
            for(int k = 0 ; k < size ; k++){
                ans[i][j] += mat_1[i][k] * mat_2[k][j];
            }
        }
    }

    return ans;

}