# Convenience library - not installed, will be linked into libmittleff
noinst_LTLIBRARIES = libquad.la
libquad_la_SOURCES = quad.c quad.h
libquad_la_CFLAGS = -fvisibility=hidden -I$(top_srcdir)/src -DBUILDING_MITTLEFF $(FLINT_CFLAGS)

# Private headers - not installed
noinst_HEADERS = quad.h
