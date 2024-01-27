void
test_sinh_0 (void)
{
    const double complex z = -1.42469226000000e+00-1.93094380000000e-01*I;
    const double complex expected = -1.92160887420070e+00-4.21900298087567e-01*I;
    const double complex computed = z*mittleff(2.0, 2.0, z*z);
    TEST_VALUE(expected, computed);
}

void
test_sinh_1 (void)
{
    const double complex z = +6.64679980000000e-01+1.46133978000000e+00*I;
    const double complex expected = +7.80741317390111e-02+1.22179750677017e+00*I;
    const double complex computed = z*mittleff(2.0, 2.0, z*z);
    TEST_VALUE(expected, computed);
}

void
test_sinh_2 (void)
{
    const double complex z = +8.94857660000000e-01+1.89827224000000e+00*I;
    const double complex expected = -3.27817271209846e-01+1.35194796777751e+00*I;
    const double complex computed = z*mittleff(2.0, 2.0, z*z);
    TEST_VALUE(expected, computed);
}

void
test_sinh_3 (void)
{
    const double complex z = -1.32190510000000e-01+2.22985315000000e+00*I;
    const double complex expected = +8.11856608972997e-02+7.97487564190490e-01*I;
    const double complex computed = z*mittleff(2.0, 2.0, z*z);
    TEST_VALUE(expected, computed);
}

void
test_sinh_4 (void)
{
    const double complex z = -1.60041154000000e+00+3.95133590000000e-01*I;
    const double complex expected = -2.19349810100817e+00+9.92523644338986e-01*I;
    const double complex computed = z*mittleff(2.0, 2.0, z*z);
    TEST_VALUE(expected, computed);
}
