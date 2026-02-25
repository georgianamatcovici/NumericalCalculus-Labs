#include <iostream>
#include <iomanip>
#include <cmath>

#define PI acos(-1)

using namespace std;

double my_tan_poly(double x)
{
    int sign = (x > 0) - (x < 0);
    x = fabs(x);

    if (x > PI / 4)
        return sign * (1.0 / my_tan_poly(PI/2 - x));

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

double my_tan_cf(double x, double eps = 1e-12)
{
    int sign = (x > 0) - (x < 0);
    x = fabs(x);

    if (x > PI / 4)
        return sign * (1.0 / my_tan_cf(PI/2 - x, eps));

    double tiny = 1e-12;

    double b0 = 1.0;
    double f = b0;
    if (fabs(f) < tiny) f = tiny;

    double C = f;
    double D = 0.0;

    double delta;
    double a, b;

    int j = 1;

    do {
        a = -x*x;
        b = 2*j + 1;

        D = b + a*D;
        if (fabs(D) < tiny) D = tiny;

        C = b + a/C;
        if (fabs(C) < tiny) C = tiny;

        D = 1.0 / D;
        delta = C * D;
        f = f * delta;

        j++;

    } while (fabs(delta - 1.0) > eps);

    return sign * x / f;
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

    // 3) aproximarea functiei tangenta
    const int N = 10000;
    double max_err_poly = 0;
    double max_err_cf = 0;

    srand(0);

    for(int i = 0; i < N; i++)
    {
        double xr = -PI/2 + (double)rand()/RAND_MAX * PI;

        double exact = tan(xr);

        double err1 = fabs(exact - my_tan_poly(xr));
        double err2 = fabs(exact - my_tan_cf(xr));

        if(err1 > max_err_poly)
            max_err_poly = err1;

        if(err2 > max_err_cf)
            max_err_cf = err2;
    }

    cout << endl << "Eroare maxima polinom: " << max_err_poly << endl;
    cout << "Eroare maxima fractii continue: " << max_err_cf << endl;

    return 0;
}
