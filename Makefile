CC = gcc
CFLAGS:= -std=c99 \
	-g \
	-Wall \
	-Wextra \
	-pedantic \
	-Wfatal-errors

VALGRIND = valgrind \
	--leak-check=full \
	--show-leak-kinds=all \
	--track-origins=yes \
	--verbose

CFLAGS_UNITY=-DUNITY_OUTPUT_COLOR -DUNITY_INCLUDE_DOUBLE -I./modules/Unity/src
CFLAGS_LOG=-DLOG_USE_COLOR -I./modules/log.c/src
LDLIBS=-lflint -lgsl -lgslcblas -lm

check: test_partition 

.PHONY: unity log test_suite

unity:
	$(CC) $(CFLAGS) $(CFLAGS_UNITY) -c ./modules/Unity/src/unity.c -o build/unity.o
log:
	$(CC) $(CFLAGS) $(CFLAGS_LOG) -c ./modules/log.c/src/log.c -o build/log.o

DEBUG = -DDEBUG $(CFLAGS_LOG)

INCDIR = -I./src/flintutils -I./src/partition -I./src

test_partition: prepare unity
	$(CC) $(CFLAGS) $(INCDIR) -c ./src/flintutils/flintutils.c -o build/flintutils.o
	$(CC) $(CFLAGS) $(INCDIR) -c ./src/partition/partition.c -o build/partition.o
	$(CC) $(CFLAGS) $(DEBUG) $(CFLAGS_UNITY) -I./src/partition -I./src/num -c tests/00-test_partition.c -o build/test_partition.o
	$(CC) $(CFLAGS) \
		build/unity.o \
		build/flintutils.o \
		build/partition.o \
		build/test_partition.o \
		-o build/test_partition.x $(LDLIBS)
	$(VALGRIND) --log-file=valgrind-test_partition.log ./build/test_partition.x

test_quad: prepare log unity
	$(CC) $(CFLAGS) -I./src -I./src/num -I./src/quad -c ./src/quad/quad.c -o build/quad.o
	$(CC) $(CFLAGS) $(DEBUG) $(CFLAGS_UNITY) -I./src/quad -I./src/num -c tests/test_quad.c -o build/test_quad.o
	$(CC) $(CFLAGS) \
		build/log.o \
		build/unity.o \
		build/quad.o \
		build/new.o \
		build/num.o \
		build/test_quad.o \
		-o build/test_quad.x $(LDLIBS)
	$(VALGRIND) --log-file=valgrind-test_quad.log ./build/test_quad.x

test_suite: prepare log unity
	$(CC) $(CFLAGS) $(INCDIR) -c ./src/flintutils/flintutils.c -o build/flintutils.o
	$(CC) $(CFLAGS) $(INCDIR) -c ./src/partition/partition.c -o build/partition.o
	$(CC) $(CFLAGS) $(INCDIR) -I./src/integrate -c ./src/integrate/integrate.c -o build/integrate.o
	$(CC) $(CFLAGS) $(INCDIR) $(DEBUG) -I./src/integrate -I./src/algorithm -c ./src/algorithm/algorithm.c -o build/algorithm.o
	$(CC) $(CFLAGS) $(INCDIR) $(DEBUG) -I./src/algorithm -I./src/mittleff -c ./src/mittleff/mittleff.c -o build/mittleff.o
	$(CC) $(CFLAGS) $(INCDIR) $(DEBUG) $(CFLAGS_UNITY) -I./src/mittleff -c tests/00-test_suite.c -o build/00-test_suite.o
	$(CC) $(CFLAGS) \
		build/log.o \
		build/unity.o \
		build/flintutils.o \
		build/algorithm.o \
		build/integrate.o \
		build/partition.o \
		build/mittleff.o \
		build/00-test_suite.o \
		-o build/test.x $(LDLIBS)
	$(VALGRIND) --log-file=valgrind-test.log ./build/test.x

.PHONY: prepare clean
clean:
	@rm -rf build && mkdir build
prepare:
	@mkdir --parents build
