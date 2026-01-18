#include <check.h>
#include <mittleff/mittleff.h>

START_TEST(test_mittleff)
{
    ck_assert_int_eq(mittleff(), 5);
}
END_TEST

int main(void)
{
  return 0;
}
