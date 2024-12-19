#include"lagrange.h"


int main(){

    srand(time(NULL));

    double  a           = 0, // концы отрезков
            b           = 6.5;
    size_t  N           = 10, // количество узлов
            M           = 100; // размер сетки для вывода
    double* x           = get_random_massive(a, b, N); // точки, по которым идет аппроксимация
    double* y           = new double[N];
    double* denominator = get_denominator(x, N);
    for(size_t i=0; i<N; ++i){
    	y[i] = func(x[i]);
    }

    ofstream data;
    data.open("data");
    data << a << " " << b << " " << N << endl;
    for(size_t i=0; i<N; ++i){
        data << x[i] << " "
             << y[i] << " ";
    }
    for(size_t i=0; i<=M; ++i){
        data << a + (b-a) * i / M << " "
             << lagrange(a + (b-a) * i / M, x, y, N, denominator) << " "
             << func(a + (b-a) * i / M) << " ";
    }
    data.close();
    system("python3 File.py");

    delete [] x;
    delete [] y;
    delete [] denominator;
    return 0;

}