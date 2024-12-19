#include"best_approx.h"

double func(double x){

	return sin(5.0*x);

}

// создает массив из M равномерно распределенных вещественных чисел от a до b
// 1-е и M-е всегда a и b
// требует отчистки памяти
double* get_uniform_mesh(double a, double b, size_t M){

	double* x = new double[M];
    for(size_t i=0; i<M; ++i){
    	x[i] = a + (b-a) * i / (M-1);
    }
    return x;

}

// создает массив из N случайных вещественных чисел от a до b
// требует отчистки памяти
double* get_random_massive(double a, double b, size_t N){

    double* x = new double[N];
    double r;
    bool flag = 0; // хранит информацию о совпадении нового случайного числа с уже сгенерированными
    for(size_t i=0; i<N; ++i){
        while(true){
            r = a + (b-a) * rand() / RAND_MAX;
            if( abs(r-a) < EPS*(b-a) ){
                    flag = 1;
                }
            for(size_t j=0; j<i; ++j){
                if( abs(r-x[j]) < EPS*(b-a) ){
                    flag = 1;
                }
            }
            if( abs(r-b) < EPS*(b-a) ){
                    flag = 1;
                }
            if(flag == 0){
                x[i] = r;
                break;
            }
            else{
                flag = 0;
            }
        }
    }
    return x;

}

double scalar(double* &x, double* &y, size_t N){

    double sum = 0;
    for(size_t i=0; i<N; ++i){
        sum += x[i]*y[i];
    }
    return sum;

}


double* get_denominator_old(double* &x, size_t N){

    double* denominator = new double[N];
    for(size_t i=0; i<N; ++i){
        double product = 1.0;
        for(size_t j=0; j<N; ++j){
            if(i==j){
                continue;
            }
            product *= x[i] - x[j];
        }
        denominator[i] = product;
    }
    return denominator;

}


double basis_func(double _x, size_t i, double* &x, size_t N, size_t K, double* &denominator){

    size_t  num     = floor( (_x - x[0]) / (x[N-1]-x[0]) ); // номер интервала, на котором _x(от нуля)
    double  result  = 0.0, 
            product = 1.0,
            x_new   = _x - (x[N-1]-x[0])*num; // локальные координаты - первый подотрезок
    bool flag       = 0;
    if(i>N-1){
        i = i + (i-1)/(N-1);
    }
    size_t bas = i%N; // индекс функции внутри интервала(от 0 до N-1)
    size_t n = i/N;
    
    if(bas == N-1){
        if( x[0] + (x[N-1]-x[0])*(n+1) < _x and _x < x[0] + (x[N-1]-x[0])*(n+2) ){
            bas = 0;
            flag = 1;
        }
        if( abs(_x - x[(n+1)*(N-1)]) < (x[K * N - K]-x[0])*EPS){return 1.0;}
    }
    if(num != n and flag == 0){ // вне нужного отрезка функции нулевые
        return 0.0;
    }

    product = 1.0;
    for(size_t j=0; j<N; ++j){
        if(bas==j){
            continue;
        }
        product *= (x_new - x[j]);
    }
    result += product / denominator[bas];
    return result;

}


// кладет в coef скалярные произведения базисных функций
//        в B скалярные произведения исходной функции на базисные
// выделяет память
// требует отчисти памяти
void get_sistem(double** &coef, size_t N, size_t K, size_t L, double* &B, double* &x, double* &denominator){

    double** mesh   = new double*[K];
    for(size_t i=0; i<K; ++i){
        mesh[i]     = get_random_massive(x[i*(N-1)], x[(i+1)*(N-1)-1], L);
    }
    double b        = 0.0;
    double sum_1    = 0.0;
    size_t M        = K * N - K + 1; // количество узлов
    coef            = new double*[M];
    B               = new double[M];
    for(size_t i=0; i<M; ++i){
        coef[i] = new double[M];
        for (size_t j=0; j<M; ++j)
        {
            coef[i][j] = 0.0;
        }
    }
    for(size_t i=0; i<M; ++i){ // перебор всех пар базисных функций
        for(size_t j=i; j<M; ++j){
            b = 0.0;
            for(size_t m=0; m<K; ++m){
                for(size_t p=0; p<L; ++p){
                    double _x = mesh[m][p];
                    sum_1 += basis_func(_x, i, x, N, K, denominator) 
                           * basis_func(_x, j, x, N, K, denominator);
                    b     += basis_func(_x, i, x, N, K, denominator)
                           * func(_x);
                }
            }
            coef[i][j] = sum_1;
            coef[j][i] = sum_1;
            sum_1 = 0.0;
        }
        B[i] = b;
    }

}

double linear_comb(double* coef, size_t L, double _x, double* &x, size_t N, size_t K, double* &denominator){

    double sum = 0;
    for(size_t i=0; i<L; ++i){
        sum += coef[i] * basis_func(_x, i, x, N, K, denominator);
    }
    return sum;

}


double l_1(double* D, size_t M){

    double result = 0;
    for(size_t i=0; i<M; ++i){
        result += abs(D[i]);
    }
    return result;

}


double l_2(double* D, size_t M){

    double result = 0;
    for(size_t i=0; i<M; ++i){
        result += pow(D[i], 2);
    }
    return pow(result, 1.0 / 2.0);

}


double l_inf(double* D, size_t M){

    double result = abs(D[0]);
    for(size_t i=1; i<M; ++i){
        if(abs(D[i]) > result)
            result = abs(D[i]);
    }
    return result;

}



// взял из интернета вместо использования библиотечной функции(метод Гаусса)
// https://teacher.ucoz.net/Lection/C/Lection6-7.pdf
int SLAU(double **matrica_a,int n,double *massiv_b, double *x)
{
    int i,j,k,r;
    double c,M,max,s, **a, *b;
    a=new double *[n];
    for(i=0;i<n;i++)
        a[i]=new double[n];
    b=new double [n];
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            a[i][j]=matrica_a[i][j];
    for(i=0;i<n;i++)
        b[i]=massiv_b[i];
    for(k=0;k<n;k++)
    {
        max=fabs(a[k][k]);
        r=k;
        for(i=k+1;i<n;i++)
        if (fabs(a[i][k])>max)
        {
            max=fabs(a[i][k]);
            r=i;
        }
        for(j=0;j<n;j++)
        {
            c=a[k][j]; a[k][j]=a[r][j];
            a[r][j]=c;
        }
        c=b[k];b[k]=b[r];b[r]=c;
        for(i=k+1;i<n;i++)
        {
            for(M=a[i][k]/a[k][k],j=k;j<n;j++)
                a[i][j]-=M*a[k][j];
        b[i]-=M*b[k];
        }
    }
    if (a[n-1][n-1]==0)
        if(b[n-1]==0)
            return -1;
        else return -2;
    else
    {
        for(i=n-1;i>=0;i--)
        {
            for(s=0,j=i+1;j<n;j++)
                s+=a[i][j]*x[j];
            x[i]=(b[i]-s)/a[i][i];
        }
        return 0;
    }
    for(i=0;i<n;i++)
    delete [] a[i];
    delete [] a;
    delete [] b;
}

