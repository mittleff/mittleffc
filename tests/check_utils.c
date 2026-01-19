/**
 * @file check_utils.c
 * @brief Unit tests for the private utils module
 */

#include <check.h>
#include "../src/utils/utils.h"

START_TEST(test_utils_add)
{
    ck_assert_int_eq(utils_add(2, 3), 5);
    ck_assert_int_eq(utils_add(0, 0), 0);
    ck_assert_int_eq(utils_add(-1, 1), 0);
    ck_assert_int_eq(utils_add(100, 200), 300);
}
END_TEST

START_TEST(test_utils_multiply)
{
    ck_assert_int_eq(utils_multiply(2, 3), 6);
    ck_assert_int_eq(utils_multiply(0, 5), 0);
    ck_assert_int_eq(utils_multiply(1, 10), 10);
    ck_assert_int_eq(utils_multiply(-2, 3), -6);
}
END_TEST

Suite *
utils_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Utils");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_utils_add);
    tcase_add_test(tc_core, test_utils_multiply);
    suite_add_tcase(s, tc_core);

    return s;
}

int
main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = utils_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? 0 : 1;
}
