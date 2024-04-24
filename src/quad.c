#include <complex.h>
#include <gsl/gsl_integration.h>

typedef double complex (*fptr)(double, void*);

double
fn_re (double x, void * params)
{
    fptr fn;
    return creal((*fn)(x));
}

double
fn_im (double x, void * params)
{
    fptr fn;
    return cimag((*fn)(x));
}

double complex
quad (double complex (*fn)(double, void *), double a, double b)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
       
    double re, im, err;
     
    gsl_function F_re, F_im;
    F_re.function = &fn_re;
    F_im.function = &fn_im;
    F_re.params = fn;
    F_im.params = fn;
     
    gsl_integration_qags (&F_re, a, b, 0, 1e-7, 1000,
                          w, &re, &error);
    gsl_integration_qags (&F_im, a, b, 0, 1e-7, 1000,
                          w, &im, &error); 
     
    /* printf ("result          = % .18f\n", result); */
    /* printf ("exact result    = % .18f\n", expected); */
    /* printf ("estimated error = % .18f\n", error); */
    /* printf ("actual error    = % .18f\n", result - expected); */
    /* printf ("intervals =  %d\n", w->size); */
     
    return re + I * im;
}

/* #include <stdio.h> */
/* #include <complex.h> */
/* #include <math.h> */

/* typedef struct { double m; double fm; double simp; } triple; */

/* /\* "structured" adaptive version, translated from Racket *\/ */
/* triple _quad_simpsons_mem(double (*f)(double), double a, double fa, double b, double fb) { */
/*     // Evaluates Simpson's Rule, also returning m and f(m) to reuse. */
/*     double m = (a + b) / 2; */
/*     double fm = f(m); */
/*     double simp = fabs(b - a) / 6 * (fa + 4*fm + fb); */
/*     triple t = {m, fm, simp}; */
/*     return t; */
/* } */

/* double _quad_asr(double (*f)(double), double a, double fa, double b, double fb, double eps, double whole, double m, double fm) { */
/*     // Efficient recursive implementation of adaptive Simpson's rule. */
/*     // Function values at the start, middle, end of the intervals are retained. */
/*     triple lt = _quad_simpsons_mem(f, a, fa, m, fm); */
/*     triple rt = _quad_simpsons_mem(f, m, fm, b, fb); */
/*     double delta = lt.simp + rt.simp - whole; */
/*     if (fabs(delta) <= eps * 15) return lt.simp + rt.simp + delta/15; */
/*     return _quad_asr(f, a, fa, m, fm, eps/2, lt.simp, lt.m, lt.fm) + */
/*            _quad_asr(f, m, fm, b, fb, eps/2, rt.simp, rt.m, rt.fm); */
/* } */

/* double quad_asr(double (*f)(double), double a, double b, double eps) { */
/*     // Integrate f from a to b using ASR with max error of eps. */
/*     double fa = f(a); */
/*     double fb = f(b); */
/*     triple t = _quad_simpsons_mem(f, a, fa, b, fb); */
/*     return _quad_asr(f, a, fa, b, fb, eps, t.simp, t.m, t.fm); */
/* } */

/* double */
/* quad (double (*fn)(double x, void * params), */
/*       double a, */
/*       double b, */
/*       double acc, */
/*       int depth, */
/*     void * params) */
/* { */
/*     // Integrate f from a to b using ASR with max error of eps. */
/*     double fa = fn(a, params); */
/*     double fb = fn(b, params); */
/*     triple t = _quad_simpsons_mem(f, a, fa, b, fb); */
/*     return _quad_asr(f, a, fa, b, fb, eps, t.simp, t.m, t.fm); */
/* } */

/* /\* int main(){ *\/ */
/* /\*     double a = 0.0, b = 1.0; *\/ */
/* /\*     double sinx = quad_asr(sin, a, b, 1e-09); *\/ */
/* /\*     printf("Simpson's integration of sine from %g to %g = %f\n", a, b, sinx); *\/ */
/* /\*     return 0; *\/ */
/* /\* } *\/ */
