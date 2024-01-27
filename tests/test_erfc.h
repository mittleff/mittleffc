void
test_erfc_0 (void)
{
    const double complex z = -7.33057219000000e-02-5.11934762000000e-01*I;
    const double complex expected = 7.25165453880409e-01-4.32914368208589e-01*I;
    const double complex computed = mittleff(0.5, 1.0, z);
    TEST_VALUE(expected, computed);
}

void
test_erfc_1 (void)
{
    const double complex z = +1.00809273000000e+01+2.22251668000000e+00*I;
    const double complex expected = 1.32220943230009e+42+1.43926327412783e+42*I;
    const double complex computed = mittleff(0.5, 1.0, z);
    TEST_VALUE(expected, computed);
}

void
test_erfc_2 (void)
{
    const double complex z = -8.81638303000000e+00+4.53794350000000e+00*I;
    const double complex expected = 5.05454404812233e-02+2.57564201381802e-02*I;
    const double complex computed = mittleff(0.5, 1.0, z);
    TEST_VALUE(expected, computed);
}

void
test_erfc_3 (void)
{
    const double complex z = -3.22342758000000e-01+8.45119872000000e+00*I;
    const double complex expected = 2.59774904698260e-03+6.71347813921331e-02*I;
    const double complex computed = mittleff(0.5, 1.0, z);
    TEST_VALUE(expected, computed);
}

void
test_erfc_4 (void)
{
    const double complex z = -3.75588680000000e-01-9.83203507000000e+00*I;
    const double complex expected = 2.22360956770567e-03-5.75980077079615e-02*I;
    const double complex computed = mittleff(0.5, 1.0, z);
    TEST_VALUE(expected, computed);
}

void
test_erfc_5 (void)
{
    const double complex z = +4.08373780000000e+00+2.53485316000000e+00*I;
    const double complex expected = -1.58178109396067e+04+5.43930514682910e+04*I;
    const double complex computed = mittleff(0.5, 1.0, z);
    TEST_VALUE(expected, computed);
}

void
test_erfc_6 (void)
{
    const double complex z = -5.00775165000000e+00+4.08876443000000e+00*I;
    const double complex expected = 6.80477811127566e-02+5.42607316062483e-02*I;
    const double complex computed = mittleff(0.5, 1.0, z);
    TEST_VALUE(expected, computed);
}


