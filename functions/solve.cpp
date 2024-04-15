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

vector<vector <double>> inverse(vector<vector <double>> &arr, const int NG) {

    vector<vector <double>> ans(NG , vector<double> (NG)); 
    alglib::real_2d_array A;
    alglib::ae_int_t info;
    alglib::matinvreport rep;

    A.setlength(NG, NG);

    double* MA = new double[NG * NG];

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

    return ans;
}

vector<double> mult(vector<vector <double>> &mat, vector<double> &arr){

    const int size = mat.size();
    int i, j;

    vector<double> ans(size,0);

    for(i = 0 ; i < size ; i++){
        for(j = 0 ; j < size ; j++){ 
            ans[i] += mat[i][j] * arr[j];
        }
    }

    return ans;
}

void mult_part(vector<vector <double>> &mat, vector<double> &arr, vector<double> &ans, const int init, const int end, const int G, const int M){
    vector<double> ans_p(M*G,0);
    int index;

    for(int g = 0 ; g < G ; g++){
        for(int i = init ; i < end ; i++){
            index = i + g*M;
            ans_p[index] = 0;
            for(int j = 0 ; j < M*G ; j++){
                ans_p[index] += mat[index][j] * arr[j];
            }
        }
    }

    for(int g = 0 ; g < G ; g++){
        for(int i = init ; i < end ; i++){
            index = i + g*M;
            ans[index] = ans_p[index];
        }
    }
}





// Function to perform matrix-vector multiplication for a subset of rows
void multiplySubset(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector,
                    std::vector<double>& result, int start, int end) {
    int cols = matrix[0].size();
    for (int i = start; i < end; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

// Function to perform matrix-vector multiplication using threads
std::vector<double> matrixVectorMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector) {
    int rows = matrix.size();
    std::vector<double> result(rows, 0.0);

    // Determine the number of threads to use (equal to number of CPU cores)
    int num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;

    // Divide the work among the threads
    int chunk_size = rows / num_threads;
    int start = 0;
    for (int i = 0; i < num_threads; ++i) {
        int end = (i == num_threads - 1) ? rows : start + chunk_size;
        threads.push_back(std::thread(multiplySubset, std::ref(matrix), std::ref(vector), std::ref(result), start, end));
        start = end;
    }

    // Join the threads
    for (std::thread& t : threads) {
        t.join();
    }

    return result;
}














void mult_part_p(vector<vector <double>> &mat, vector<double> &arr, vector<double> &ans, const int init, const int end, const int G, const int M){
    vector<double> ans_p(M*G,0);
    int index;

    for(int g = 0 ; g < G ; g++){
        for(int i = init ; i < end ; i++){
            index = i + g*M;
            ans_p[index] = 0;
            #pragma omp parallel for 
            for(int j = 0 ; j < M*G ; j++){
                ans_p[index] += mat[index][j] * arr[j];
            }
        }
    }

    for(int g = 0 ; g < G ; g++){
        for(int i = init ; i < end ; i++){
            index = i + g*M;
            ans[index] = ans_p[index];
        }
    }
}

vector<double> mult_sub(vector<vector <double>> &mat, vector<double> &arr){

    const int size = mat.size();

    vector<double> ans(size,0);

    #pragma omp parallel 
    {
        vector<double> ans_p(size,0);
        #pragma omp for collapse(2)
        for(int i = 0 ; i < size ; i++){
            for(int j = 0 ; j < size ; j++){
                ans_p[i] += mat[i][j] * arr[j];
            }
        }
        #pragma omp critical
        {
            for(int i = 0 ; i < size ; i++) ans[i] -= ans_p[i];
        }
    
    }
    
    return ans;
}

vector<vector<double>> mult_mat(vector<vector <double>> &mat_1, vector<vector<double>> &mat_2){

    const int size = mat_1.size();
    vector<vector<double>> ans (size, vector<double> (size));

    for(int i = 0 ; i < size ; i++){
        for(int k = 0 ; k < size ; k++){
            for(int j = 0 ; j < size ; j++){
                ans[i][j] += mat_1[i][k] * mat_2[k][j];
            }
        }
    }

    return ans;

}

vector<vector<double>> mult_mat_sub(vector<vector <double>> &mat_1, vector<vector<double>> &mat_2){

    const int size = mat_1.size();
    vector<vector<double>> ans (size, vector<double> (size));

    for(int i = 0 ; i < size ; i++){
        for(int j = 0 ; j < size ; j++){
            for(int k = 0 ; k < size ; k++){
                ans[i][j] -= mat_1[i][k] * mat_2[k][j];
            }
        }
    }

    return ans;

}

vector<double> sum_arr(vector<double> &arr_1, vector<double> &arr_2){

    const int size = arr_1.size();
    vector<double> ans(size,0);

    for(int i = 0 ; i < size ; i++){
        ans[i] = arr_1[i] + arr_2[i];
    }

    return ans;
}

vector<double> parallelAddition (unsigned N, vector<double> &A, vector<double> &B) {
    unsigned i;

    vector<double> C(N, 0);

    // #pragma omp critical    
    for (i = 0; i < N; ++i) C[i] = A[i] + B[i];
        
    return C;
}

vector<double> parallelSubtraction (unsigned N, vector<double> &A, vector<double> &B) {
    unsigned i;

    vector<double> C(N, 0);

    // #pragma omp critical    
    for (i = 0; i < N; ++i) C[i] = A[i] - B[i];
        
    return C;
}

vector<double> sub_arr(vector<double> &arr_1, vector<double> &arr_2){

    const int size = arr_1.size();
    vector<double> ans(size,0);

    
        for(int i = 0 ; i < size ; i++){
            ans[i] = arr_1[i] - arr_2[i];
        }    
    

    return ans;
}

vector<vector<double>> sum_mat(vector<vector<double>> &mat_1, vector<vector<double>> &mat_2){

    const int size = mat_1[0].size();
    vector<vector <double>> ans(size , vector<double> (size));

    for(int i = 0 ; i < size ; i++){
        for(int j = 0 ; j < size ; j++){
            ans[i][j] = mat_1[i][j] + mat_2[i][j];   
        }
    }    

    return ans;
}

vector<vector<double>> sub_mat(vector<vector<double>> &mat_1, vector<vector<double>> &mat_2){

    const int size = mat_1[0].size();
    vector<vector <double>> ans(size , vector<double> (size));

    for(int i = 0 ; i < size ; i++){
        for(int j = 0 ; j < size ; j++){
            ans[i][j] = mat_1[i][j] - mat_2[i][j];   
        }
    }    

    return ans;
}

vector<vector<double>> add_eye(vector<vector<double>> &mat_1){

    const int size = mat_1[0].size();
    vector<vector <double>> ans(size , vector<double> (size));

    for(int i = 0 ; i < size ; i++){
        for(int j = 0 ; j < size ; j++){
            ans[i][j] = mat_1[i][j];
        }
    }

    for(int i = 0 ; i < size ; i++){
        ans[i][i] = mat_1[i][i] - 1;   
    }    

    return ans;
}


