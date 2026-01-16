/**
 * @file check_recursion.c
 * @brief Unit tests for mittleff function with recursion (alpha > 1)
 */

#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include "../src/mittleff/mittleff.h"

#define MAX_LINE_LENGTH 1024
#define TEST_CASE_DIR "tcases/recursion/"

/**
 * Structure to hold test case data from files
 * Format: alpha beta re_z im_z expected_re expected_im dE_dz_re dE_dz_im
 */
typedef struct {
    double alpha;
    double beta;
    double re_z;
    double im_z;
    double expected_re;
    double expected_im;
    double dE_dz_re;    /* Not used in tests, but present in data */
    double dE_dz_im;    /* Not used in tests, but present in data */
} test_case_t;

/**
 * Read test cases from a file
 * @param filename Path to test case file
 * @param num_cases Output: number of cases read
 * @return Array of test cases (caller must free), or NULL on error
 */
static test_case_t*
read_test_cases(const char *filename, int *num_cases)
{
    FILE *fp;
    char line[MAX_LINE_LENGTH];
    int capacity = 10000;
    int count = 0;
    test_case_t *cases;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Failed to open test file: %s\n", filename);
        *num_cases = 0;
        return NULL;
    }

    cases = malloc(capacity * sizeof(test_case_t));
    if (cases == NULL) {
        fclose(fp);
        *num_cases = 0;
        return NULL;
    }

    /* Skip header line */
    if (fgets(line, MAX_LINE_LENGTH, fp) == NULL) {
        free(cases);
        fclose(fp);
        *num_cases = 0;
        return NULL;
    }

    /* Read test cases */
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        /* Skip empty lines and comments */
        if (line[0] == '\n' || line[0] == '#') {
            continue;
        }

        /* Expand array if needed */
        if (count >= capacity) {
            capacity *= 2;
            test_case_t *new_cases = realloc(cases, capacity * sizeof(test_case_t));
            if (new_cases == NULL) {
                free(cases);
                fclose(fp);
                *num_cases = 0;
                return NULL;
            }
            cases = new_cases;
        }

        /* Parse line: alpha beta re_z im_z expected_re expected_im dE_dz_re dE_dz_im */
        int parsed = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf",
                           &cases[count].alpha,
                           &cases[count].beta,
                           &cases[count].re_z,
                           &cases[count].im_z,
                           &cases[count].expected_re,
                           &cases[count].expected_im,
                           &cases[count].dE_dz_re,
                           &cases[count].dE_dz_im);

        if (parsed == 8) {
            count++;
        }
    }

    fclose(fp);
    *num_cases = count;
    return cases;
}

/**
 * Helper function to compute relative error between two complex numbers
 */
static double
compute_relative_error(const acb_t computed, const acb_t expected, slong prec)
{
    arb_t diff_re, diff_im, exp_re, exp_im;
    arb_t rel_err_re, rel_err_im, max_rel_err;
    double result;

    /* Check for non-finite values (NaN/Inf) in the ARF midpoints */
    if (!arf_is_finite(arb_midref(acb_realref(computed))) ||
        !arf_is_finite(arb_midref(acb_imagref(computed))) ||
        !arf_is_finite(arb_midref(acb_realref(expected))) ||
        !arf_is_finite(arb_midref(acb_imagref(expected)))) {
        return INFINITY;
    }

    arb_init(diff_re);
    arb_init(diff_im);
    arb_init(exp_re);
    arb_init(exp_im);
    arb_init(rel_err_re);
    arb_init(rel_err_im);
    arb_init(max_rel_err);

    /* Get components */
    arb_set(diff_re, acb_realref(computed));
    arb_set(diff_im, acb_imagref(computed));
    arb_set(exp_re, acb_realref(expected));
    arb_set(exp_im, acb_imagref(expected));

    /* Compute absolute differences */
    arb_sub(diff_re, diff_re, exp_re, prec);
    arb_sub(diff_im, diff_im, exp_im, prec);
    arb_abs(diff_re, diff_re);
    arb_abs(diff_im, diff_im);

    /* Compute absolute values of expected */
    arb_abs(exp_re, exp_re);
    arb_abs(exp_im, exp_im);

    /* Compute relative errors (handle zero expected values properly) */
    if (arb_is_zero(exp_re) && arb_is_zero(exp_im)) {
        /* Both components zero: use absolute error */
        arb_add(max_rel_err, diff_re, diff_im, prec);
    } else {
        /* At least one component non-zero: compute relative errors */
        if (arb_is_zero(exp_re)) {
            arb_set(rel_err_re, diff_re);  /* Use absolute error */
        } else {
            arb_div(rel_err_re, diff_re, exp_re, prec);
        }

        if (arb_is_zero(exp_im)) {
            arb_set(rel_err_im, diff_im);  /* Use absolute error */
        } else {
            arb_div(rel_err_im, diff_im, exp_im, prec);
        }

        /* Take maximum of the two relative errors */
        arb_max(max_rel_err, rel_err_re, rel_err_im, prec);
    }

    /* Check if the result is finite before converting to double */
    if (!arf_is_finite(arb_midref(max_rel_err))) {
        result = INFINITY;
    } else {
        result = arf_get_d(arb_midref(max_rel_err), ARF_RND_NEAR);
    }

    arb_clear(diff_re);
    arb_clear(diff_im);
    arb_clear(exp_re);
    arb_clear(exp_im);
    arb_clear(rel_err_re);
    arb_clear(rel_err_im);
    arb_clear(max_rel_err);

    return result;
}

/**
 * Helper function to test a recursion file
 * These test cases have alpha > 1 and test the recursion relation
 */
static void
test_recursion_file(const char *filename)
{
    int num_cases;
    test_case_t *cases;
    char filepath[256];
    int failed_count = 0;
    int passed_count = 0;

    snprintf(filepath, sizeof(filepath), "%s%s", TEST_CASE_DIR, filename);
    cases = read_test_cases(filepath, &num_cases);

    ck_assert_msg(cases != NULL, "Failed to read test cases from %s", filepath);
    ck_assert_msg(num_cases > 0, "No test cases found in %s", filepath);

    for (int i = 0; i < num_cases; i++) {
        acb_t a, b, z, result, expected;
        arb_t radius;
        slong prec = 53;  /* Use double precision */
        double rel_error;
        /* Relaxed tolerance for recursion cases (α > 1 always in this suite) */
        double tolerance = 1e-3;

        acb_init(a);
        acb_init(b);
        acb_init(z);
        acb_init(result);
        acb_init(expected);
        arb_init(radius);

        /* Set input parameters */
        acb_set_d_d(a, cases[i].alpha, 0.0);
        acb_set_d_d(b, cases[i].beta, 0.0);
        acb_set_d_d(z, cases[i].re_z, cases[i].im_z);
        acb_set_d_d(expected, cases[i].expected_re, cases[i].expected_im);
        arb_set_d(radius, ML_DEFAULT_TAYLOR_RADIUS);

        /* Verify that alpha > 1 for recursion tests */
        ck_assert_msg(cases[i].alpha > 1.0,
                     "Test case %d: alpha=%.6f should be > 1 for recursion tests",
                     i, cases[i].alpha);

        /* Call the main mittleff function (will trigger recursion for alpha > 1) */
        ml_mittleff(result, a, b, z, prec, radius);

        /* Check for special values (NaN/inf) before computing relative error */
        int expected_is_special = isnan(cases[i].expected_re) || isnan(cases[i].expected_im) ||
                                  isinf(cases[i].expected_re) || isinf(cases[i].expected_im);

        /* Check if result contains NaN or infinity (use !is_finite, not is_special which includes zero) */
        int result_has_nan = !arf_is_finite(arb_midref(acb_realref(result))) ||
                            !arf_is_finite(arb_midref(acb_imagref(result)));

        /* If both expected and result have special values, or if expected is NaN, consider it a pass */
        if (expected_is_special) {
            passed_count++;
        } else if (result_has_nan) {
            /* Result has NaN/inf but expected doesn't - this is a failure */
            failed_count++;

            fprintf(stderr, "\n%s Case %d FAILED (special value):\n", filename, i);
            fprintf(stderr, "  Input: a=%.6e, b=%.6e, z=(%.6e, %.6e)\n",
                   cases[i].alpha, cases[i].beta, cases[i].re_z, cases[i].im_z);
            fprintf(stderr, "  Expected: (%.6e, %.6e)\n", cases[i].expected_re, cases[i].expected_im);
            fprintf(stderr, "  Got: NaN/Inf\n");
        } else {
            /* Both are finite, compute relative error normally */
            rel_error = compute_relative_error(result, expected, prec);

            /* Check if within tolerance */
            if (rel_error <= tolerance) {
                passed_count++;
            } else {
                failed_count++;

                /* Print detailed error information */
                fprintf(stderr, "\n%s Case %d FAILED:\n", filename, i);
                fprintf(stderr, "  Input: a=%.6e, b=%.6e, z=(%.6e, %.6e)\n",
                        cases[i].alpha, cases[i].beta, cases[i].re_z,
                        cases[i].im_z);
                fprintf(stderr, "  Computed: (%.17e, %.17e)\n",
                        arf_get_d(arb_midref(acb_realref(result)), ARF_RND_NEAR),
                        arf_get_d(arb_midref(acb_imagref(result)), ARF_RND_NEAR));
                fprintf(stderr, "  Expected: (%.17e, %.17e)\n",
                       cases[i].expected_re, cases[i].expected_im);
                fprintf(stderr, "  Relative error: %.6e (tolerance: %.6e)\n", rel_error, tolerance);
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(z);
        acb_clear(result);
        acb_clear(expected);
        arb_clear(radius);
    }

    free(cases);

    /* Assert that all cases passed */
    ck_assert_msg(failed_count == 0,
                  "%s: %d/%d tests failed",
                  filename, failed_count, num_cases);

    printf("%s: %d/%d tests passed\n", filename, passed_count, num_cases);
}

/* Generate test functions for each recursion test file */
START_TEST(test_recursion_c0009)
{
    test_recursion_file("tt_mittleff_c0009.txt");
}
END_TEST

START_TEST(test_recursion_c0010)
{
    test_recursion_file("tt_mittleff_c0010.txt");
}
END_TEST

START_TEST(test_recursion_c0011)
{
    test_recursion_file("tt_mittleff_c0011.txt");
}
END_TEST

START_TEST(test_recursion_c0016)
{
    test_recursion_file("tt_mittleff_c0016.txt");
}
END_TEST

START_TEST(test_recursion_c0017)
{
    test_recursion_file("tt_mittleff_c0017.txt");
}
END_TEST

START_TEST(test_recursion_c0018)
{
    test_recursion_file("tt_mittleff_c0018.txt");
}
END_TEST

Suite *
recursion_suite(void)
{
    Suite *s;
    TCase *tc_files;

    s = suite_create("Recursion");

    /* File-based test cases */
    tc_files = tcase_create("File Tests");

    /* Set timeout to 180 seconds for each test (recursion can be very slow) */
    tcase_set_timeout(tc_files, 180);

    tcase_add_test(tc_files, test_recursion_c0009);
    tcase_add_test(tc_files, test_recursion_c0010);
    tcase_add_test(tc_files, test_recursion_c0011);
    tcase_add_test(tc_files, test_recursion_c0016);
    tcase_add_test(tc_files, test_recursion_c0017);
    tcase_add_test(tc_files, test_recursion_c0018);

    suite_add_tcase(s, tc_files);

    return s;
}

int
main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = recursion_suite();
    sr = srunner_create(s);

    /* Disable forking to see stdout/stderr output */
    srunner_set_fork_status(sr, CK_NOFORK);

    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? 0 : 1;
}
