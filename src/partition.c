#include <complex.h>
#include <math.h>

int open_wedgep (double complex z,
                 double phi1,
                 double phi2);
int closed_wedgep (double complex z,
                 double phi1,
                 double phi2);

double
compute_r1 (double a,
    double acc)
{
    double c0;
    
    c0 = (pow(1.3, 1 - a))/(M_PI * sin(M_PI * a));
    return pow(-2 * log(acc / c0), a);
}

int
in_region_G1 (double complex z,
              double a)
{
    double delta , phi1, phi2;

    delta = 0.125 * M_PI * a;
    phi1 = - M_PI * a + delta;
    phi2 = + M_PI * a - delta;

    return open_wedgep(z, phi1, phi2);
}

int
in_region_G2 (double complex z,
              double a)
{
    double delta, tau, phi1, phi2;

    delta = 0.125 * M_PI * a;
    tau = 0.5 * M_PI * (a + 1);
    tau = (delta < tau) ? delta : tau;

    phi1 = + M_PI * a + tau;
    phi2 = - phi1;

    return open_wedgep(z, phi1, phi2);
}

int
in_region_G3 (double complex z,
              double a)
{
    double delta, tau, phi1, phi2;

    delta = 0.125 * M_PI * a;
    tau = 0.5 * M_PI * (a + 1);
    tau = (delta < tau) ? delta : tau;

    phi1 = + M_PI * a - delta;
    phi2 = + M_PI * a + tau;

    return closed_wedgep(z, phi1, phi2);
}

int
in_region_G4 (double complex z,
              double a)
{
    double delta, tau, phi1, phi2;

    delta = 0.125 * M_PI * a;
    tau = 0.5 * M_PI * (a + 1);
    tau = (delta < tau) ? delta : tau;

    phi1 = - M_PI * a - tau;
    phi2 = - M_PI * a + delta;

    return closed_wedgep(z, phi1, phi2);
}

int
in_region_G5 (double complex z,
              double a,
              double r0)
{
    double phi1, phi2;
    
    phi1 = - (5.0/6.0) * M_PI * a;
    phi2 = + (5.0/6.0) * M_PI * a;

    return closed_wedgep(z, phi1, phi2) && (cabs(z) >= r0);
}

int
in_region_G6 (double complex z,
              double a,
              double r0)
{
    double phi1, phi2;
    
    phi1 = + (5.0/6.0) * M_PI * a;
    phi2 = - (5.0/6.0) * M_PI * a;

    return open_wedgep(z, phi1, phi2) && (cabs(z) >= r0);
}

double
cross (double x, double y)
{
    return sin(x - y);
}

int
clockwisep (double a, double b, int closed)
{
    if (closed)
        return cross(a, b) >= 0;
    else
        return cross(a, b) > 0;
}

int between (double complex z,
             double phi1,
             double phi2,
             int closed)
{
    double th;

    th = carg(z);

    if (clockwisep(phi1, phi2, closed))
        return !between(z, phi2, phi1, closed);
    else
        return ((!clockwisep(phi1, th, closed)) &&
                (clockwisep(phi2, th, closed)));
}

int open_wedgep (double complex z,
                 double phi1,
                 double phi2)
{
    return between(z, phi1, phi2, 0);
}

int closed_wedgep (double complex z,
                 double phi1,
                 double phi2)
{
    return between(z, phi1, phi2, 1);
}
              
