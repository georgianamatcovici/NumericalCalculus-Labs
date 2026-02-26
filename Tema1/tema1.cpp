#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>

#define PI acos(-1)

using namespace std;

double reduce_x(double x) {
    x = fmod(x, PI);
    if (x > PI / 2.0) x -= PI;
    if (x < -PI / 2.0) x += PI;
    return x;
}

double my_tan_poly(double x)
{
    x = reduce_x(x);
    if (fabs(fabs(x) - PI/2.0) < 1e-15) return INFINITY;
    int sign = (x >= 0) - (x < 0);
    x = fabs(x);

    if (x > PI / 4)
        return sign * (1.0 / my_tan_poly(PI/2.0 - x));

    double c1 = 1.0/3.0;
    double c2 = 2.0/15.0;
    double c3 = 17.0/315.0;
    double c4 = 62.0/2835.0;

    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x2*x2;
    double x6 = x4*x2;

    double result = x + x3*(c1 + c2*x2 + c3*x4 + c4*x6);

    return sign * result;
}

double my_tan_cf(double x, double eps)
{
    x = reduce_x(x);
    if (fabs(fabs(x) - PI/2.0) < 1e-15) return INFINITY;
    if (fabs(x) < 1e-15) return 0.0;

    double tiny = 1e-12;

    double f = tiny;
    double C = f;
    double D = 0.0;
    double delta;

    for (int j = 1; j < 1000; j++) {
        double a = (j == 1) ? x : -x*x;
        double b = 2.0 * j - 1.0;

        D = b + a * D;
        if (fabs(D) < 1e-18) D = tiny;
        D = 1.0 / D;

        C = b + a / C;
        if (fabs(C) < 1e-18) C = tiny;

        delta = C * D;
        f *= delta;

        if (fabs(delta - 1.0) < eps) break;
    }

    return f;
}

int main() {
    cout << fixed << setprecision(20);

    // 1) Determinare precizie masina

    double u = 1.0;
    int m = 0;

    while (1.0 + u != 1.0) {
        u /= 10.0;
        m++;
    }

    u *= 10.0;
    m -= 1;

    cout << "Precizia masina:" << endl;
    cout << "u = 10^(-" << m << ") = " << u << endl;
    cout << "1 + u      = " << (1.0 + u) << endl;
    cout << "1 + u/10   = " << (1.0 + u/10.0) << endl;
    cout << "((1+u) - 1) = " << ((1.0 + u) - 1.0) << endl;
    cout << "((1+u/10) - 1) = " << ((1.0 + u/10.0) - 1.0) << endl;

    // 2) Neasociativitate adunare
    double x = 1.0;
    double y = u / 10.0;
    double z = u / 10.0;

    double t1 = x + y;
    double r1 = t1 + z;

    double t2 = y + z;
    double r2 = x + t2;

    cout << endl << "Adunare:" << endl;
    cout << "(x+y)+z = " << r1 << endl;
    cout << "x+(y+z) = " << r2 << endl;

    if (r1 != r2)
        cout << "Adunarea NU este asociativa." << endl;
    else
        cout << "Nu se observa diferenta in acest caz." << endl;

    // 2) Neasociativitate inmultire

    x = 1e308;
    y = 1e308;
    z = 1e-308;

    double t3 = x * y;
    double a  = t3 * z;

    double t4 = y * z;
    double b  = x * t4;

    cout << endl << "Inmultire:" << endl;
    cout << "(x*y)*z = " << a << endl;
    cout << "x*(y*z) = " << b << endl;

    if (a != b)
        cout << "Inmultirea NU este asociativa." << endl;
    else
        cout << "Nu se observa diferenta in acest caz." << endl;

    // 3) Aproximarea functiei tangenta
    const int N = 10000;
    double x_vals[N];
    srand(time(0));

    for(int i = 0; i < N; i++)
        x_vals[i] = -PI/2.0 + (double)rand()/RAND_MAX * PI;

    int p;
    cout << "\nIntrodu p pentru precizie (eps = 10^-p): ";
    cin >> p;
    double eps = pow(10.0, -p);

    // masurare polinom
    auto start_poly = chrono::high_resolution_clock::now();
    double max_err_poly = 0;
    for(int i = 0; i < N; i++) {
        double exact = tan(x_vals[i]);
        double approx = my_tan_poly(x_vals[i]);
        if(isfinite(exact) && isfinite(approx))
            max_err_poly = max(max_err_poly, fabs(exact - approx));
    }
    auto end_poly = chrono::high_resolution_clock::now();

    // masurare fractii continue
    auto start_cf = chrono::high_resolution_clock::now();
    double max_err_cf = 0;
    for(int i = 0; i < N; i++) {
        double exact = tan(x_vals[i]);
        double approx = my_tan_cf(x_vals[i], eps);
        if(isfinite(exact) && isfinite(approx))
            max_err_cf = max(max_err_cf, fabs(exact - approx));
    }
    auto end_cf = chrono::high_resolution_clock::now();

    chrono::duration<double, milli> time_poly = end_poly - start_poly;
    chrono::duration<double, milli> time_cf = end_cf - start_cf;

    cout << fixed << setprecision(15);
    cout << endl << "REZULTATE 10.000 ITERATII" << endl;
    cout << "Polinom: Max Err = " << max_err_poly << " | Timp = " << time_poly.count() << " ms" << endl;
    cout << "Fractii: Max Err = " << max_err_cf << " | Timp = " << time_cf.count() << " ms" << endl;
    return 0;
}
