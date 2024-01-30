void
test_siam_1 (void)
{
    const double complex computed = mittleff(0.6, 0.8, 7.0);
    const double complex expected = +4.24680224735076e+11;
    TEST_VALUE(expected, computed);
}


void
test_siam_2 (void)
{
    const double complex computed = mittleff(0.6, 0.8, 20.0);
    const double complex expected = +4.50513132816147e+64;        
    TEST_VALUE(expected, computed);
}

void
test_siam_3 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex computed = mittleff(0.6, 0.8, -7.0);
    const double complex expected = 0.0364029650876388;
    TEST_VALUE(expected, computed);
}

void
test_siam_4 (void)
{
    const double complex computed = mittleff(0.6, 0.8, -50.0);
    const double complex expected = 0.0044638678216643;
    TEST_VALUE(expected, computed);
}

void
test_siam_5 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex computed = mittleff(0.6, 0.8, -2.16311896062463+6.65739561406608*I);
    const double complex expected = 0.00509750816218177+0.0329981074690976*I;    
    TEST_VALUE(expected, computed);
}

void
test_siam_6 (void)
{
    const double complex computed = mittleff(0.6, 0.8, -6.18033988749895+19.0211303259031*I);
    const double complex expected = 0.00282134530403973+0.0107554765459201*I;
    TEST_VALUE(expected, computed);
}

void
test_siam_7 (void)
{
    const double complex computed = mittleff(0.6, 1.25, 7.0);
    const double complex expected = 98682128538.563;
    TEST_VALUE(expected, computed);
}

void
test_siam_8 (void)
{        
    const double complex computed = mittleff(0.6, 1.25, 20.0);
    const double complex expected = 4.76359640442376e+63+9.21339224649432e-19*I;
    TEST_VALUE(expected, computed);
}

void
test_siam_9 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex computed = mittleff(0.6, 1.25, -7.0);
    const double complex expected = 0.101261033685572;
    TEST_VALUE(expected, computed);
}

void
test_siam_10 (void)
{
    const double complex computed = mittleff(0.6, 1.25, -50.0);
    const double complex expected = 0.0144197663438114-7.6778268720786e-20*I;
    TEST_VALUE(expected, computed);
}

void
test_siam_11 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex computed = mittleff(0.6, 1.25, -2.16311896062463+6.65739561406608*I);
    const double complex expected = 0.0333902562082633+0.0980431639835736*I;
    TEST_VALUE(expected, computed);
}

void
test_siam_12 (void)
{
   
    const double complex expected = 0.011289456355613+0.0342852434746551*I;
    const double complex computed = mittleff(0.6, 1.25,-6.18033988749895+19.0211303259031*I);
    TEST_VALUE(expected, computed);
}

void
test_siam_13 (void)
{
    const double complex expected = 76147703794042.9;
    const double complex computed = mittleff(0.6, -0.8, 7.0);
    TEST_VALUE(expected, computed);
}

void
test_siam_14 (void)
{
    const double complex expected = 1.32776365747668e+68;
    const double complex computed = mittleff(0.6, -0.8, 20.0);
    TEST_VALUE(expected, computed);
}

void
test_siam_15 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex expected = 0.0501291913319244;
    const double complex computed = mittleff(0.6, -0.8, -7.0);
    TEST_VALUE(expected, computed);
}

void
test_siam_16 (void)
{
    const double complex expected = 0.00751163297262774;
    const double complex computed = mittleff(0.6, -0.8, -50);
    TEST_VALUE(expected, computed);
}

void
test_siam_17 (void)
{
    //TEST_IGNORE_MESSAGE("Check integration routine (region G6)");
    const double complex expected = 0.0193182614473201+0.0537209282676945*I;
    const double complex computed = mittleff(0.6, -0.8, -2.16311896062463+6.65739561406608*I);
    TEST_VALUE(expected, computed);
}

void
test_siam_18 (void)
{
    const double complex expected = 0.00592228306634142+0.0179734030934208*I;
    const double complex computed = mittleff(0.6, -0.8, -6.18033988749895+19.0211303259031*I);
    TEST_VALUE(expected, computed);
}
