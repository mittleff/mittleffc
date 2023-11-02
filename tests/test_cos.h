void
test_cos_0 (void)
{
    const double complex z = -1.42469226000000e+00-1.93094380000000e-01*I;
    const double complex expected = +1.48307362593645e-01-1.92226474311045e-01*I;
    const double complex computed = mittleff(2.0, 1.0, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_cos_1 (void)
{
    const double complex z = +6.64679980000000e-01+1.46133978000000e+00*I;
    const double complex expected = +1.78818881219233e+00-1.25822734251227e+00*I;
    const double complex computed = mittleff(2.0, 1.0, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_cos_2 (void)
{
    const double complex z = +8.94857660000000e-01+1.89827224000000e+00*I;
    const double complex expected = +2.13470542364750e+00-2.54495655003334e+00*I;
    const double complex computed = mittleff(2.0, 1.0, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_cos_3 (void)
{
    const double complex z = -1.32190510000000e-01+2.22985315000000e+00*I;
    const double complex expected = +4.66199107535284e+00+6.05710930000310e-01*I;
    const double complex computed = mittleff(2.0, 1.0, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_cos_4 (void)
{
    const double complex z = -1.60041154000000e+00+3.95133590000000e-01*I;
    const double complex expected = -3.19526988999422e-02+4.05318417916536e-01*I;
    const double complex computed = mittleff(2.0, 1.0, -z*z);
    TEST_VALUE(expected, computed);
}
