CC:=gcc
SHELL:=/bin/bash
CFLAGS=-Wall -Werror
LDFLAGS=-lm -lflint

.PHONY: clean

shared: build/libmittleff.so

SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:.c=.o)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< 

build/libmittleff.so: $(OBJS)
	$(CC) -fPIC -shared $(CFLAGS) $^ -o $@ $(LDFLAGS)

build/libpartition.so: src/partition.o
	$(CC) -fPIC -shared $(CFLAGS) $^ -o $@ $(LDFLAGS)
build/libalgorithm.so: src/rgamma.o src/partition.o src/algorithm.o
	$(CC) -fPIC -shared $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm src/*.o

###########
# Testing #
###########
.PHONY: test test-partition test-algorithm test-mittleff

test: test-partition test-algorithm test-mittleff test-recursion

test-partition: build/libpartition.so
	time LD_LIBRARY_PATH=build/ GUILE_LOAD_PATH=./ guile -s tests/test-partition.scm

test-algorithm: build/libalgorithm.so
	time LD_LIBRARY_PATH=build/ GUILE_LOAD_PATH=./ guile -s tests/test-algorithm.scm

test-mittleff: shared
	time LD_LIBRARY_PATH=build/ GUILE_LOAD_PATH=./ guile -s tests/test-mittleff.scm

test-recursion: shared
	time LD_LIBRARY_PATH=build/ GUILE_LOAD_PATH=./ guile -s tests/test-recursion.scm
