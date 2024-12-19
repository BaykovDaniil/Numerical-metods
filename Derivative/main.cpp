#include"derivative.h"


int main(){
    double a            = -6.0, // концы отрезков
           b            = 5.0;
    size_t M            = 20, // количество узлов
           Mviz         = 100;
    double h            = (b-a) / (M-1); // шаг сетки 1
    double* x1          = get_uniform_mesh(a, b, M);
    double* x2          = get_uniform_mesh(a, b, 2*M-1);
    double* D1          = Derivative(x1, M);
    double* D2          = Derivative(x2, 2*M-1);
    double* D3          = Runge(D1, D2, M);

    double* D1_error    = new double[M];
    double* D2_error    = new double[2*M-1];
    double* D3_error    = new double[M];
    double* R_error     = TheoreticalRungeError(D1, D2, M);

    for(size_t i=0; i<M; ++i)
        D1_error[i] = D1[i] - D_func(x1[i]);
    for(size_t i=0; i<2*M-1; ++i)
        D2_error[i] = D2[i] - D_func(x2[i]);
    for(size_t i=0; i<M; ++i)
        D3_error[i] = D3[i] - D_func(x1[i]);

    cout << setprecision(6);
    cout << scientific;
    cout << endl;
    cout << "Таблица 1 - абсолютная погрешность"
         << endl;
    cout << setw(10) << "l_1"
         << setw(14) << "l_2"
         << setw(16) << "l_inf" << endl;
    cout << "h"
         << setw(18) << l_1(D1_error, M)
         << setw(14) << l_2(D1_error, M)
         << setw(14) << l_inf(D1_error, M) << endl;
    cout << "h/2"
         << setw(16) << l_1(D2_error, 2*M-1)
         << setw(14) << l_2(D2_error, 2*M-1)
         << setw(14) << l_inf(D2_error, 2*M-1) << endl;
    cout << "Runge"
         << setw(14) << l_1(D3_error, M)
         << setw(14) << l_2(D3_error, M)
         << setw(14) << l_inf(D3_error, M) << endl;
    cout << endl;

    cout << "Таблица 2 - относительная погрешность"
         << endl;
    cout << setw(10) << "l_1"
         << setw(14) << "l_2"
         << setw(16) << "l_inf" << endl;
    cout << "h"
         << setw(18) << l_1(D1_error, M) / l_1(D1, M)
         << setw(14) << l_2(D1_error, M) / l_2(D1, M)
         << setw(14) << l_inf(D1_error, M) / l_inf(D1, M) << endl;
    cout << "h/2"
         << setw(16) << l_1(D2_error, 2*M-1) / l_1(D2, 2*M-1)
         << setw(14) << l_2(D2_error, 2*M-1) / l_2(D2, 2*M-1)
         << setw(14) << l_inf(D2_error, 2*M-1) / l_inf(D2, 2*M-1) << endl;
    cout << "Runge"
         << setw(14) << l_1(D3_error, M) / l_1(D3, M)
         << setw(14) << l_2(D3_error, M) / l_2(D3, M)
         << setw(14) << l_inf(D3_error, M) / l_inf(D3, M) << endl;
    cout << endl;

    cout << "Таблица 3 - погрешность Рунге" 
         << endl;
    cout << setw(16) << "l_1"
         << setw(14) << "l_2"
         << setw(16) << "l_inf" << endl;
    cout << "Real_error"
         << setw(15) << l_1(D1_error, M)
         << setw(14) << l_2(D1_error, M)
         << setw(14) << l_inf(D1_error, M) << endl;
    cout << "Runge_error"
         << setw(14) << l_1(R_error, M)
         << setw(14) << l_2(R_error, M)
         << setw(14) << l_inf(R_error, M) << endl;
    cout << endl;

    ofstream data;
    data.open("data");
    data << a << " " << b << " " << M << endl;
    for(size_t i=0; i<M; ++i){
        data << x1[i] << " " 
             << D1[i] << " ";
    }
    data << endl;
    for(size_t i=0; i<2*M-1; ++i){
        data << x2[i] << " "
             << D2[i] << " ";
    }
    data << endl;
    for(size_t i=0; i<M; ++i){
        data << x1[i] << " " 
             << D3[i] << " ";
    }
    data << endl;
    for(size_t i=0; i<=Mviz; ++i){
        data << a + (b-a) * i / Mviz << " "
             << D_func(a + (b-a) * i / Mviz) << " ";
    }
    data.close();
    system("python3 File.py");


    delete [] x1;
    delete [] x2;
    delete [] D1;
    delete [] D2;
    delete [] D3;
    delete [] D1_error;
    delete [] D2_error;
    delete [] D3_error;
    delete [] R_error;

    return 0;
}