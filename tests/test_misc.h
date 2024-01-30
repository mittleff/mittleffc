void
test_misc_1 (void)
{
    const double complex computed = mittleff(0.5, 0.5, 0.5);
    const double complex expected = 1.5403698281390346;
    TEST_VALUE(expected, computed);
}

void
test_misc_2 (void)
{
    const double complex computed = mittleff(1.5, 0.5, 0.5);
    const double complex expected = 1.1448466286155243;
    TEST_VALUE(expected, computed);
}

void
test_misc_3 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex computed = mittleff(2.3, 1.0, 0.7+2.0*I);
    const double complex expected = 1.201890136368392+0.7895394560075035*I;
    TEST_VALUE(expected, computed);
}

void
test_misc_4 (void)
{
    const double complex computed = mittleff(2.3, 1.0, 0.7+0.2*I);
    const double complex expected = 1.268233154873853+0.07914994421659409*I;
    TEST_VALUE(expected, computed);
}
void
test_misc_5 (void)
{
    TEST_IGNORE_MESSAGE("Not handling bignum yet");
    const double complex computed = mittleff(0.3, 1.0, 100.0);
    const double complex expected = 8.721285946907744692995882256235296113802695745418015206361825134909144332670706e+2015816;
    TEST_VALUE(expected, computed);
}

void
test_misc_6 (void)
{
    const double complex computed = mittleff(0.9, 0.5, 22.0+22.0*I);
    const double complex expected = -2.7808021618204008e13-2.8561425165239754e13*I;
    TEST_VALUE(expected, computed);
}

void
test_misc_7 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex computed = mittleff(0.1, 1.05, 0.9+0.5*I);
    const double complex expected = 0.17617901349590603+2.063981943021305*I;
    TEST_VALUE(expected, computed);
}

void
test_misc_8 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex computed = mittleff(4.1, 1.0, 1.0);
    const double complex expected = 1.0358176744122032;
    TEST_VALUE(expected, computed);
}

void
test_misc_9 (void)
{
    const double complex computed = mittleff(0.5, 1.0, -12.0);
    const double complex expected = 0.046854221014893775;
    TEST_VALUE(expected, computed);
}

void
test_misc_10 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex computed = mittleff(0.125, 1.0, -1.0);
    const double complex expected = 0.481952081535048487353320281623;
    TEST_VALUE(expected, computed);
}
