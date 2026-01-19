# Convenience library - not installed, will be linked into libmittleff
noinst_LTLIBRARIES = libalgorithm.la
libalgorithm_la_SOURCES = algorithm.c algorithm.h
libalgorithm_la_CFLAGS = -fvisibility=hidden -I$(top_srcdir)/src -DBUILDING_MITTLEFF $(FLINT_CFLAGS)

# Private headers - not installed
noinst_HEADERS = algorithm.h
