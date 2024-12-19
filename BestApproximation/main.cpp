#include"best_approx.h"


int main(){
    srand(time(NULL));
    // концы отрезков
    double a            = 1.0;
    double b            = 3.0;
    size_t K            = 4; // количество интервалов
    size_t N            = 4; // количество точек на интервале
    size_t L            = 100; // количество случайных точек
    size_t M            = K * N - K + 1; // количество узлов
    double h            = (b-a) / (M-1); // шаг сетки
    size_t Mviz         = 100;
    double* x           = get_uniform_mesh(a, b, M);
    double* denominator = get_denominator_old(x, N);
    double** coef       = nullptr;
    double* B           = nullptr;
    double* ans_of_sis  = new double[K*N];
    
    get_sistem(coef, N, K, L, B, x, denominator);
    int res = SLAU(coef, M, B, ans_of_sis);

    double* error = new double[100*M];
    double* approx = new double[100*M];
    for(size_t i=0; i<100*M; ++i){
        approx[i] = linear_comb(ans_of_sis, M, a + (b-a) * i / (100*M-1), x, N, K, denominator);
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
             << func(x[i]) << " ";
    }
    data << endl;

    for(size_t i=0; i<=Mviz; ++i){
        data << a + (b-a) * i / Mviz << " "
             << linear_comb(ans_of_sis, M, a + (b-a) * i / Mviz, x, N, K, denominator) << " ";
    }
    data << endl;

    for(size_t i=0; i<=Mviz; ++i){
        data << a + (b-a) * i / Mviz << " "
             << func(a + (b-a) * i / Mviz) << " ";
    }

    data.close();
    system("python3 File.py");

    delete [] x;
    delete [] denominator;
    for(size_t i=0; i<M; ++i){
        delete [] coef[i];
    }
    delete [] coef;
    delete [] B;
    delete [] ans_of_sis;
    delete [] error;
    delete [] approx;
    return 0;
}