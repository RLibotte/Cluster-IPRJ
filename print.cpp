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

void printArray(std::vector<double> arr){
    for(int i = 0 ; i < arr.size() ; i++){
        printf("%2.6e ", arr[i]);
    }
    printf("\n\n");
}

void print2dArray(std::vector<std::vector<double>> arr){
    for(int i = 0 ; i < arr.size() ; i++){
        for(int j = 0 ; j < arr[0].size() ; j++){
            printf("%2.6e\t", arr[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void print3dArray(std::vector<std::vector<std::vector<double>>> arr){
    for(int i = 0 ; i < arr.size() ; i++){
        for(int j = 0 ; j < arr[0].size() ; j++){
            for(int k = 0 ; k < arr[0][0].size() ; k++){
                printf("%2.6e\t", arr[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n\n");
}