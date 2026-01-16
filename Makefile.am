SUBDIRS = src . tests examples
dist_doc_DATA = README AUTHORS NEWS ChangeLog COPYING

# Install pkg-config file
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mittleff.pc

# Custom cleaning rules
clean-local:
	git clean -fX
