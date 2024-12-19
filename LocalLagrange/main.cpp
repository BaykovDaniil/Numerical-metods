#include"lagrange.h"


int main(){
    // концы отрезков
    double a            = 1.0;
    double b            = 3.0;
    size_t K            = 5; // количество интервалов
    size_t N            = 4; // степень многочленов
    size_t M            = K * N - K + 1; // количество узлов интерполяции
    double h            = (b-a) / (M-1); // шаг сетки
    size_t Mviz         = 100;
    double* x           = get_uniform_mesh(a, b, M);
    double* y           = new double[M];
    double* denominator = get_denominator(x, N);
    for(int i=0; i<M; ++i){
        y[i] = func(x[i]);
    }

    double* error = new double[100*M];
    double* approx = new double[100*M];
    for(size_t i=0; i<100*M; ++i){
        approx[i] = lagrange(a + (b-a) * i / (100*M-1), x, y, N, denominator);
        error[i]  = func(a + (b-a) * i / (100*M-1)) - 
                    approx[i];
    }

    cout<<setprecision(6);
    cout << scientific << endl;
    cout << "Таблица 1 - погрешность по разным нормам" 
         << endl;
    cout << setw(19) << "l_1"
         << setw(14) << "l_2"
         << setw(16) << "l_inf" << endl;
    cout << "Absolute_error"
         << setw(14) << l_1(error, 100*M)
         << setw(14) << l_2(error, 100*M)
         << setw(14) << l_inf(error, 100*M) << endl;
    cout << "Relative_error"
         << setw(14) << l_1(error, 100*M) / l_1(approx, 100*M)
         << setw(14) << l_2(error, 100*M) / l_2(approx, 100*M)
         << setw(14) << l_inf(error, 100*M) / l_inf(approx, 100*M) << endl;
    cout << endl;

    ofstream data;
    data.open("data");
    data << a << " " << b << " " << M << " " << N << endl;
    for(size_t i=0; i<M; ++i){
        data << x[i] << " "
             << y[i] << " ";
    }
    for(size_t i=0; i<=Mviz; ++i){
        data << a + (b-a) * i / Mviz << " "
             << lagrange(a + (b-a) * i / Mviz, x, y, N, denominator) << " "
             << func(a + (b-a) * i / Mviz) << " ";
    }
    data.close();
    system("python3 File.py");

    delete [] x;
    delete [] y;
    delete [] denominator;
    delete [] error;
    delete [] approx;

    return 0;
}