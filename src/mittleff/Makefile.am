lib_LTLIBRARIES = libmittleff.la
libmittleff_la_SOURCES = mittleff.c
libmittleff_la_LDFLAGS = -version-info 0:0:0
libmittleff_la_CFLAGS = -fvisibility=hidden -I$(top_srcdir)/src -DBUILDING_MITTLEFF $(FLINT_CFLAGS)
libmittleff_la_LIBADD = $(top_builddir)/src/utils/libutils.la \
                        $(top_builddir)/src/partition/libpartition.la \
                        $(top_builddir)/src/algorithm/libalgorithm.la \
                        $(FLINT_LIBS) -lflint

# Install headers to $(includedir)/mittleff/
mittleffincludedir = $(includedir)/mittleff
mittleffinclude_HEADERS = mittleff.h
