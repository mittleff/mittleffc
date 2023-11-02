

void
test_sin_0 (void)
{
    const double complex z = -1.42469226000000e+00-1.93094380000000e-01*I;
    const double complex expected = -1.00784724885197e+00-2.82866292071741e-02*I;
    const double complex computed = z*mittleff(2.0, 2.0, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_sin_1 (void)
{
    const double complex z = +6.64679980000000e-01+1.46133978000000e+00*I;
    const double complex expected = +1.40128057728241e+00+1.60563708193143e+00*I;
    const double complex computed = z*mittleff(2.0, 2.0, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_sin_2 (void)
{
    const double complex z = +8.94857660000000e-01+1.89827224000000e+00*I;
    const double complex expected = +2.66183979862789e+00+2.04096901440268e+00*I;
    const double complex computed = z*mittleff(2.0, 2.0, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_sin_3 (void)
{
    const double complex z = -1.32190510000000e-01+2.22985315000000e+00*I;
    const double complex expected = -6.19885888230589e-01+4.55538511767942e+00*I;
    const double complex computed = z*mittleff(2.0, 2.0, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_sin_4 (void)
{
    const double complex z = -1.60041154000000e+00+3.95133590000000e-01*I;
    const double complex expected = -1.07861309811772e+00-1.20071018874968e-02*I;
    const double complex computed = z*mittleff(2.0, 2.0, -z*z);
    TEST_VALUE(expected, computed);
}
