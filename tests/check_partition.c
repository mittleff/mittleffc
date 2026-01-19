/**
 * @file check_partition.c
 * @brief Unit tests for the partition module
 */

#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include "../src/partition/partition.h"

#define MAX_LINE_LENGTH 1024
#define TEST_CASE_DIR "tcases/partition/"

/**
 * Structure to hold test case data from files
 */
typedef struct {
    double r0;
    double alpha;
    double re_z;
    double im_z;
    int expected_region;
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

        /* Parse line: r0 alpha re_z im_z region */
        int parsed = sscanf(line, "%lf %lf %lf %lf %d",
                           &cases[count].r0,
                           &cases[count].alpha,
                           &cases[count].re_z,
                           &cases[count].im_z,
                           &cases[count].expected_region);

        if (parsed == 5) {
            count++;
        }
    }

    fclose(fp);
    *num_cases = count;
    return cases;
}

/**
 * Helper function to test a partition file
 */
static void
test_partition_file(const char *filename)
{
  int num_cases;
  test_case_t *cases;
  char filepath[256];

  snprintf(filepath, sizeof(filepath), "%s%s", TEST_CASE_DIR, filename);
  cases = read_test_cases(filepath, &num_cases);

  ck_assert_msg(cases != NULL, "Failed to read test cases from %s", filepath);
  ck_assert_msg(num_cases > 0, "No test cases found in %s", filepath);

  for (int i = 0; i < num_cases; i++) {
    acb_t alpha;
    acb_t z;
    arb_t taylor_radius;
    slong prec = 48;  /* Match Scheme test generation precision */

    acb_init(z);
    acb_init(alpha);
    arb_init(taylor_radius);
        
    acb_set_d_d(z, cases[i].re_z, cases[i].im_z);
    acb_set_d_d(alpha, cases[i].alpha, 0);
    arb_set_d(taylor_radius, cases[i].r0);             

    /* Call the partition function */
    int result = ml_compute_region(z, alpha, taylor_radius, prec);

    /* Compare with expected result */
    ck_assert_msg(result == cases[i].expected_region,
                  "%s Case %d failed: expected %d, got %d (r0=%.6e, alpha=%.6e, re_z=%.6e, im_z=%.6e)",
                  filename, i, cases[i].expected_region, result,
                  cases[i].r0, cases[i].alpha, cases[i].re_z, cases[i].im_z);

    acb_clear(z);
    acb_clear(alpha);
    arb_clear(taylor_radius);
  }

  free(cases);
}

START_TEST(test_partition_c0001)
{
    test_partition_file("tt_partition_c0001.txt");
}
END_TEST

START_TEST(test_partition_c0002)
{
    test_partition_file("tt_partition_c0002.txt");
}
END_TEST

START_TEST(test_partition_c0003)
{
    test_partition_file("tt_partition_c0003.txt");
}
END_TEST

START_TEST(test_partition_c0004)
{
    test_partition_file("tt_partition_c0004.txt");
}
END_TEST

START_TEST(test_partition_c0005)
{
    test_partition_file("tt_partition_c0005.txt");
}
END_TEST

START_TEST(test_partition_c0006)
{
    test_partition_file("tt_partition_c0006.txt");
}
END_TEST

START_TEST(test_partition_c0007)
{
    test_partition_file("tt_partition_c0007.txt");
}
END_TEST

START_TEST(test_partition_c0008)
{
    test_partition_file("tt_partition_c0008.txt");
}
END_TEST

START_TEST(test_partition_c0009)
{
    test_partition_file("tt_partition_c0009.txt");
}
END_TEST

START_TEST(test_partition_c0010)
{
    test_partition_file("tt_partition_c0010.txt");
}
END_TEST

START_TEST(test_partition_c0011)
{
    test_partition_file("tt_partition_c0011.txt");
}
END_TEST

START_TEST(test_partition_c0012)
{
    test_partition_file("tt_partition_c0012.txt");
}
END_TEST

/**
 * Test basic ml_clockwise functionality
 */
/* START_TEST(test_ml_clockwise_basic) */
/* { */
/*     arb_t a, b; */
/*     arb_init(a); */
/*     arb_init(b); */

/*     /\* Test: 90 degrees is clockwise from 0 degrees *\/ */
/*     arb_set_d(a, 1.5708);  /\* ~π/2 radians *\/ */
/*     arb_set_d(b, 0.0); */
/*     ck_assert_int_eq(ml_clockwise(a, b, 0), 1); */

/*     /\* Test: 0 degrees is not clockwise from 90 degrees *\/ */
/*     arb_set_d(a, 0.0); */
/*     arb_set_d(b, 1.5708); */
/*     ck_assert_int_eq(ml_clockwise(a, b, 0), 0); */

/*     arb_clear(a); */
/*     arb_clear(b); */
/* } */
/* END_TEST */

Suite *
partition_suite(void)
{
    Suite *s;
    /* TCase *tc_core; */
    TCase *tc_files;

    s = suite_create("Partition");

    /* Core functionality tests */
    /* tc_core = tcase_create("Core"); */
    /* tcase_add_test(tc_core, test_ml_clockwise_basic); */
    /* suite_add_tcase(s, tc_core); */

    /* File-based test cases */
    tc_files = tcase_create("File Tests");
    tcase_add_test(tc_files, test_partition_c0001);
    tcase_add_test(tc_files, test_partition_c0002);
    tcase_add_test(tc_files, test_partition_c0003);
    tcase_add_test(tc_files, test_partition_c0004);
    tcase_add_test(tc_files, test_partition_c0005);
    tcase_add_test(tc_files, test_partition_c0006);
    tcase_add_test(tc_files, test_partition_c0007);
    tcase_add_test(tc_files, test_partition_c0008);
    tcase_add_test(tc_files, test_partition_c0009);
    tcase_add_test(tc_files, test_partition_c0010);
    tcase_add_test(tc_files, test_partition_c0011);
    tcase_add_test(tc_files, test_partition_c0012);
    /* Add more test files as needed:
     * tcase_add_test(tc_files, test_partition_c0003);
     * tcase_add_test(tc_files, test_partition_c0004);
     * ...
     */
    suite_add_tcase(s, tc_files);

    return s;
}

int
main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = partition_suite();
    sr = srunner_create(s);

    /* Disable forking to see stdout/stderr output */
    srunner_set_fork_status(sr, CK_NOFORK);

    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? 0 : 1;
}
