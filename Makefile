TARGET ?= libmittleff.so

CC = gcc
CFLAGS:= -std=c99 \
	-Wall \
	-Wextra \
	-pedantic \
	-Wfatal-errors

# For logging
# LOGGING_CFLAGS = -DLOG_USE_COLOR  -DLOGLEVEL=$(LOGLEVEL)
# LOGGING_SRCS=./logging/src/log.c
# LOGGING_INCDIR=./logging/src

# For Unity (testing library)
UNITY_CFLAGS = -DUNITY_INCLUDE_DOUBLE
UNITY_SRCS=./modules/unity/src/unity.c
UNITY_INCDIR=./modules/unity/src/

# For numeric
# NUMERIC_CFLAGS = -DM_PI=3.14159265358979323846 -D_TOLERANCE=1e-7
# NUMERIC_SRCS=$(shell find ./numeric/src/ -type f -name '*.c')
# NUMERIC_INCDIR=./numeric/src/
# NUMERIC_LDFLAGS = -lm

# For integration
# INTEGRATION_CFLAGS =
# INTEGRATION_SRCS=$(shell find ./integration/src/ -type f -name '*.c')
# INTEGRATION_INCDIR=./integration/src/

# For sf
# SF_CFLAGS =
# SF_SRCS=$(shell find ./sf/src/ -type f -name '*.c')
# SF_INCDIR=./sf/src/
# SF_LDFLAGS = -lgsl -lgslcblas -larb

# For mittleff
MITTLEFF_CFLAGS = -fPIC -O2 -g
MITTLEFF_SRCS=$(shell find ./src/ -type f -name '*.c')
MITTLEFF_INCDIR=./include/

# INCDIR = $(LOGGING_INCDIR) $(UNITY_INCDIR) $(NUMERIC_INCDIR) $(INTEGRATION_INCDIR) $(SF_INCDIR) $(MITTLEFF_INCDIR)
# INCFLAGS=$(foreach d,$(INCDIR),-I$d)
# CFLAGS+=$(LOGGING_CFLAGS) $(UNITY_CFLAGS) $(NUMERIC_CFLAGS) $(INTEGRATION_CFLAGS) $(SF_CFLAGS) $(MITTLEFF_CFLAGS)
# LDFLAGS+=$(NUMERIC_LDFLAGS) $(SF_LDFLAGS)

INCDIR = $(UNITY_INCDIR) $(MITTLEFF_INCDIR)
INCFLAGS=$(foreach d,$(INCDIR),-I$d)
CFLAGS+=$(UNITY_CFLAGS) $(MITTLEFF_CFLAGS)

.PHONY: test
SRCS_TEST := $(UNITY_SRCS) $(MITTLEFF_SRCS) ./test/test.c
OBJS_TEST := $(SRCS_TEST:%.c=%.o)

test: test.out
	./test.out

test.out: $(OBJS_TEST)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(INCFLAGS) $(CFLAGS) $(CCFLAGS) -c $^ -o $@

all: test


	#LDFLAGS = -shared
#LDMATH = -lm -lgsl -lgslcblas -larb
#LOGFLAGS = -DLOG_USE_COLOR
#INCFLAGS = -I./include
#SRCS := $(shell find src/ -type f -name '*.c')
#OBJS := $(SRCS:%.c=%.o)
#TESTS_DIR:=tests
#
#%.o: %.c
#	$(CC) $(INCFLAGS) $(CFLAGS) $(LOGFLAGS) -c $^ -o $@
#
#.PHONY: all
#all: $(TARGET)
#
#$(TARGET): $(OBJS)
#	$(CC) $^ -o $@ $(LDFLAGS) $(LDMATH)
#
#.PHONY: doc
#doc: Doxyfile
#	mkdir -p doc/
#	doxygen
#
.PHONY: clean clean-all
clean:
	find . -iname "*.o" -type f -delete
#clean-all: clean
#	find . -iname "libmittleff.so" -type f -delete
#	rm -rf tests/tests_num.c
