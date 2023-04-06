////////////////////////////////////////////////////////////////
//         Universidade do Estado do Rio de Janeiro           //
//                           PPGMC                            //
//                                                            //
//              Simulador Transporte de Nêutrons              //
//                                                            // 
//               Autor: Rafael Barbosa Libotte                //
////////////////////////////////////////////////////////////////

using namespace std;

#define PI 3.14159265358979323846  

vector<double> insertion_sort(int N, vector<double> array){
    double max, min, size;
    for(int i = 1 ; i < N ; i++){
        for (int j = 0; j < i; j++) {
            if(array[i] < array[j]){
                min = array[i];
                max = array[j];
                array[j] = min;
                array[i] = max;
            }
        }
    }
    return array;
}

void get_GQ(int N, vector<double> &mi, vector<double> &w){

    mi.resize(N);
    w.resize(N);

    double M = 0.5*(N+1);
    double u,u1,P1,P2,P3,DP;

    for (int i = 0 ; i < M-1 ; i++){
        u = cos(PI * ((i+1) - 0.25)/(N+0.5));
            do {
                P1 = 1;
                P2 = 0;
                for (int j = 0 ; j < N; j++){
                    P3 = P2;
                    P2 = P1;
                    P1 = ((2 * j + 1) * u * P2 - j * P3)/(j+1);
                }
                DP = N*(u*P1 - P2)/(pow(u, 2)-1);
                u1 = u;
                u = u1-(P1/DP);
            } while (fabs(u-u1) > pow(10, -15));
            
        mi[i] = u;
        mi[N-i-1] = -u;
        w[N-i-1] = 2/((1 - pow(u,2)) * pow(DP,2));
        w[i] = w[N-i-1];
    }

    mi = insertion_sort(N, mi);

     for(int i = 0; i < N/2 ; i++){
        mi[i] = mi[i+N/2];
        mi[i+N/2] = -mi[i+N/2];
    }
    
    for(int i = 0; i < N/2 ; i++){
        w[i] = w[i+N/2];
    }
}
