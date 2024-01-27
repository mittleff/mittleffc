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

check: test_suite

.PHONY: num unity log test_suite
num:
	$(CC) $(CFLAGS) -I./src -I./src/num -c ./src/num/new.c -o build/new.o
	$(CC) $(CFLAGS) -I./src -I./src/num -c ./src/num/num.c -o build/num.o
unity:
	$(CC) $(CFLAGS) $(CFLAGS_UNITY) -c ./modules/Unity/src/unity.c -o build/unity.o
log:
	$(CC) $(CFLAGS) $(CFLAGS_LOG) -c ./modules/log.c/src/log.c -o build/log.o

DEBUG = -DDEBUG $(CFLAGS_LOG)

test_suite: prepare log num unity
	$(CC) $(CFLAGS) -I./src -I./src/num -I./src/partition -c ./src/partition/partition.c -o build/partition.o
	$(CC) $(CFLAGS) -I./src -I./src/num -I./src/quad -c ./src/quad/quad.c -o build/quad.o
	$(CC) $(CFLAGS) -I./src -I./src/num -I./src/quad -I./src/integrate -c ./src/integrate/integrate.c -o build/integrate.o
	$(CC) $(CFLAGS) $(DEBUG) -I./src -I./src/num -I./src/integrate -I./src/algorithm -c ./src/algorithm/algorithm.c -o build/algorithm.o
	$(CC) $(CFLAGS) $(DEBUG) -I./src/num -I./src/partition -I./src/algorithm -I./src/mittleff -c ./src/mittleff/mittleff.c -o build/mittleff.o
	$(CC) $(CFLAGS) $(DEBUG) $(CFLAGS_UNITY) -I./src/num -I./src/partition -I./src/mittleff -c tests/00-test_suite.c -o build/00-test_suite.o
	$(CC) $(CFLAGS) \
		build/log.o \
		build/unity.o \
		build/algorithm.o \
		build/quad.o \
		build/integrate.o \
		build/new.o \
		build/num.o \
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
