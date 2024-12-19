#include"integrate.h"


int main(){
    double a    = 1.0; // концы отрезков
    double b    = 5.0;
    size_t K    = 20; // количество интервалов разбиения
    size_t M    = K+1;
    double h    = (b-a) / K; // шаг сетки
    double* x   = get_uniform_mesh(a, b, M);

    double* Int_trap        = Int_Trap(x, M);
    double* Int_rect        = Int_Rect(x, M);
    double* Int_simpson     = Int_Simpson(x, M);
    double* Int_newton      = Int_Newton(x, M);
    double* Int_gauss       = Int_Gauss(x, M);
    double True_integral    = True_Integral(a, b);

    cout << endl << scientific;
    cout << "Истина" 
         << setw(21) 
         << True_integral << endl;
    cout << setw(27) << "сетка h" 
         << setw(20) << "сетка h/2" << endl;
    cout << "Трапеции" 
         << setw(19) 
         << Int_trap[0] << ' ' 
         << Int_trap[1] << endl;
    cout << "Прямоугольники" 
         << setw(13) 
         << Int_rect[0] << ' ' 
         << Int_rect[1] << endl;
    cout << "Симпсон" 
         << setw(20) 
         << Int_simpson[0] << ' ' 
         << Int_simpson[1] << endl;
    cout << "Ньютон-Котес" 
         << setw(15) 
         << Int_newton[0] << ' ' 
         << Int_newton[1] << endl;
    cout << "Гаусс" 
         << setw(22) 
         << Int_gauss[0] << ' ' 
         << Int_gauss[1] << endl;

    cout << endl << "Таблица относительных погрешностей(%)" << endl;
    cout << setw(27) << "сетка h" 
         << setw(20) << "сетка h/2" << endl;
    cout << "Трапеции" 
         << setw(19) 
         << relative_error(True_integral, Int_trap[0]) << ' ' 
         << relative_error(True_integral, Int_trap[1]) << endl;
    cout << "Прямоугольники" 
         << setw(13) 
         << relative_error(True_integral, Int_rect[0]) << ' ' 
         << relative_error(True_integral, Int_rect[1]) << endl;
    cout << "Симпсон"
         << setw(20) 
         << relative_error(True_integral, Int_simpson[0]) << ' ' 
         << relative_error(True_integral, Int_simpson[1]) << endl;
    cout << "Ньютон-Котес"
         << setw(15) 
         << relative_error(True_integral, Int_newton[0]) << ' ' 
         << relative_error(True_integral, Int_newton[1]) << endl;
    cout << "Гаусс"
         << setw(22) 
         << relative_error(True_integral, Int_gauss[0]) << ' ' 
         << relative_error(True_integral, Int_gauss[1]) << endl;
    cout << endl;
    
    delete [] x;
    delete [] Int_trap;
    delete [] Int_rect;
    delete [] Int_simpson;
    delete [] Int_newton;
    delete [] Int_gauss;

    return 0;
}