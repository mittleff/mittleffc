CC = gcc
CFLAGS:= -std=c99 \
	-Wall \
	-Wextra \
	-pedantic \
	-Wfatal-errors

VALGRIND = valgrind \
	--leak-check=full \
	--show-leak-kinds=all \
	--track-origins=yes \
	--verbose

all: test_num

.PHONY: test_num
test_num:
	@mkdir --parents build
	$(CC) $(CFLAGS) -I./modules/Unity/src -c ./modules/Unity/src/unity.c -o build/unity.o
	$(CC) $(CFLAGS) -I./src/num -c ./src/num/new.c -o build/new.o
	$(CC) $(CFLAGS) -I./src/num -c ./src/num/num.c -o build/num.o
	$(CC) $(CFLAGS) -I./modules/Unity/src -I./src/num -DUNITY_INCLUDE_DOUBLE -c test/test_num.c -o build/test_num.o
	$(CC) $(CFLAGS) \
		build/unity.o \
		build/new.o \
		build/num.o \
		build/test_num.o \
		-o build/test_num.x -lflint -lgsl -lgslcblas -lm
	$(VALGRIND) --log-file=valgrind-test_num.txt ./build/test_num.x

DEBUG = -DDEBUG -DLOG_USE_COLOR -I./modules/log.c/src
test_mittleff:
	@mkdir --parents build
	$(CC) $(CFLAGS) -DLOG_USE_COLOR -I./modules/log.c/src         -c ./modules/log.c/src/log.c -o build/log.o
	$(CC) $(CFLAGS) -D UNITY_OUTPUT_COLOR -I./modules/Unity/src         -c ./modules/Unity/src/unity.c -DUNITY_INCLUDE_DOUBLE -o build/unity.o
	$(CC) $(CFLAGS) -I./src/num                   -c ./src/num/new.c             -o build/new.o
	$(CC) $(CFLAGS) -I./src/num                   -c ./src/num/num.c             -o build/num.o
	$(CC) $(CFLAGS) -I./src/num -I./src/partition -c ./src/partition/partition.c -o build/partition.o
	$(CC) $(CFLAGS) -I./src -I./src/num -I./src/quad      -c ./src/quad/quad.c           -o build/quad.o
	$(CC) $(CFLAGS) -I./src -I./src/num -I./src/quad -I./src/integrate -c ./src/integrate/integrate.c -o build/integrate.o
	$(CC) $(CFLAGS) $(DEBUG) -I./src -I./src/num -I./src/integrate -I./src/algorithm -c ./src/algorithm/algorithm.c -o build/algorithm.o
	$(CC) $(CFLAGS) $(DEBUG) -I./src/num -I./src/partition -I./src/algorithm -I./src/mittleff  -c ./src/mittleff/mittleff.c   -o build/mittleff.o
	$(CC) $(CFLAGS) $(DEBUG) -DDEBUG -DLOG_USE_COLOR -D UNITY_OUTPUT_COLOR -I./modules/log.c/src -I./modules/Unity/src -I./src/num -I./src/mittleff -DUNITY_INCLUDE_DOUBLE -c tests/test_mittleff.c -o build/test_mittleff.o
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
		build/test_mittleff.o \
		-o build/test_mittleff.x -lflint -lgsl -lgslcblas -lm
	$(VALGRIND) --log-file=valgrind-test_mittleff.txt ./build/test_mittleff.x


.PHONY: python_package
python_module:
	$(CC) $(CFLAGS) -fPIC -I./src/num -c ./src/num/new.c -o build/new.o
	$(CC) $(CFLAGS) -fPIC -I./src/num -c ./src/num/num.c -o build/num.o
	$(CC) $(CFLAGS) -fPIC -I./src/num -I./src/partition -c ./src/partition/partition.c -o build/partition.o
	$(CC) $(CFLAGS) -fPIC -I./src -I./src/num -I./src/quad -c ./src/quad/quad.c -o build/quad.o
	$(CC) $(CFLAGS) -fPIC -I./src -I./src/num -I./src/quad -I./src/integrate -c ./src/integrate/integrate.c -o build/integrate.o
	$(CC) $(CFLAGS) -fPIC -I./src -I./src/num -I./src/integrate -I./src/algorithm -c ./src/algorithm/algorithm.c -o build/algorithm.o
	$(CC) $(CFLAGS) -fPIC -I./src/num -I./src/partition -I./src/algorithm -I./src/mittleff  -c ./src/mittleff/mittleff.c -o build/mittleff.o
	$(CC) $(CFLAGS) -fPIC -fdiagnostics-color=always -MD -MQ -D_FILE_OFFSET_BITS=64 -Winvalid-pch -O3 -I/usr/include/python3.11  -I./src/num -I./src/mittleff -c python/mittleff_module.c -o build/mittleff_module.o
	$(CC) $(CFLAGS) -fPIC -shared \
		build/algorithm.o \
		build/quad.o \
		build/integrate.o \
		build/new.o \
		build/num.o \
		build/partition.o \
		build/mittleff.o \
		build/mittleff_module.o \
		-o build/mittleff.so -lflint -lgsl -lgslcblas -lm

# DEV_CFLAGS:=-std=c99
# DEV_INCDIR:=-I./include -I./src/include

LOGLEVEL ?= 1
all: test
.PHONY: test
test:
	@mkdir --parents build
	$(CC) $(CFLAGS) -DLOG_USE_COLOR -I./modules/log.c/src         -c ./modules/log.c/src/log.c -o build/log.o
	$(CC) $(CFLAGS) -D UNITY_OUTPUT_COLOR -I./modules/Unity/src         -c ./modules/Unity/src/unity.c -DUNITY_INCLUDE_DOUBLE -o build/unity.o
	$(CC) $(CFLAGS) -I./src/num                   -c ./src/num/new.c             -o build/new.o
	$(CC) $(CFLAGS) -I./src/num                   -c ./src/num/num.c             -o build/num.o
	$(CC) $(CFLAGS) -I./src/num -I./src/partition -c ./src/partition/partition.c -o build/partition.o
	$(CC) $(CFLAGS) -I./src -I./src/num -I./src/quad      -c ./src/quad/quad.c           -o build/quad.o
	$(CC) $(CFLAGS) -I./src -I./src/num -I./src/quad -I./src/integrate -c ./src/integrate/integrate.c -o build/integrate.o
	$(CC) $(CFLAGS) $(DEBUG) -I./src -I./src/num -I./src/integrate -I./src/algorithm -c ./src/algorithm/algorithm.c -o build/algorithm.o
	$(CC) $(CFLAGS) $(DEBUG) -I./src/num -I./src/partition -I./src/algorithm -I./src/mittleff  -c ./src/mittleff/mittleff.c   -o build/mittleff.o
	$(CC) $(CFLAGS) $(DEBUG) -DDEBUG -DLOG_USE_COLOR -D UNITY_OUTPUT_COLOR -I./modules/log.c/src -I./modules/Unity/src -I./src/num -I./src/num -I./src/partition -I./src/mittleff -DUNITY_INCLUDE_DOUBLE -c tests/00-test_suite.c -o build/00-test_suite.o
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
		-o build/test_suite.x -lflint -lgsl -lgslcblas -lm
	$(VALGRIND)	--log-file=valgrind-out.txt ./build/test_suite.x



# .PHONY: clean clean-all
# clean:
# 	@rm -rf build && mkdir build

