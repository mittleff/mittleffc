# Examples are not installed, just built for demonstration
noinst_PROGRAMS = simple

simple_SOURCES = simple.c
simple_CFLAGS = -I$(top_srcdir)/src
simple_LDADD = $(top_builddir)/src/mittleff/libmittleff.la

# Include examples in distribution
EXTRA_DIST = README.md
