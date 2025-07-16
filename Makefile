SHELL := /bin/bash
CC=gcc
CFLAGS_PROD=-Wall -Werror -fPIC -shared -DDEBUG
CFLAGS=-fPIC -shared -DDEBUG 
LDFLAGS=-lm -lflint

shared: build/lib/libutils.so build/lib/libpartition.so

build/lib/libutils.so: src/utils.c
	mkdir --parents build/{obj,lib}
	$(CC) $(CFLAGS) -I./src/include -c src/utils.c -o build/obj/utils.o
	$(CC) $(CFLAGS) build/obj/utils.o -o build/lib/libutils.so $(LDFLAGS)

build/lib/libpartition.so: src/partition.c src/utils.c
	mkdir --parents build/{obj,lib}
	$(CC) $(CFLAGS) -I./src/include -c src/utils.c -o build/obj/utils.o
	$(CC) $(CFLAGS) -I./src/include -c src/partition.c -o build/obj/partition.o
	$(CC) $(CFLAGS) build/obj/utils.o build/obj/partition.o -o build/lib/libpartition.so $(LDFLAGS)

###################
# Executable builds #
###################
.PHONY: mittleff mittleff-flint run run-flint clean-executables partition

mittleff: bin/mittleff

bin/mittleff: src/mittleff.c src/partition.c
	mkdir -p bin
	$(CC) -Wall -O2 -o $@ $^ $(LDFLAGS)

run: bin/mittleff
	./bin/mittleff

libpartition: build/lib/libpartition.so
libutils: build/lib/libutils.so

clean-executables:
	rm -f bin/mittleff

###########
# Testing #
###########
.PHONY: test test-partition test-algorithm test-mittleff

# test: test-partition test-algorithm test-mittleff test-recursion
test: test-partition

test-utils: libutils
	time LD_LIBRARY_PATH=build/lib/ GUILE_LOAD_PATH=./bindings/guile:./ guile -s tests/test-utils.scm

test-partition: libpartition
	time LD_LIBRARY_PATH=build/lib/ GUILE_LOAD_PATH=./bindings/guile:./ guile -s tests/test-partition.scm

test-algorithm: shared
	time LD_LIBRARY_PATH=build/ GUILE_LOAD_PATH=./bindings/guile:./ guile -s tests/test-algorithm.scm

test-mittleff: shared
	time LD_LIBRARY_PATH=build/ GUILE_LOAD_PATH=./bindings/guile:./ guile -s tests/test-mittleff.scm

test-recursion: shared
	time LD_LIBRARY_PATH=build/ GUILE_LOAD_PATH=./bindings/guile:./ guile -s tests/test-recursion.scm

repl:
	LD_LIBRARY_PATH=build/lib GUILE_LOAD_PATH=./bindings/guile:./ guile 

##########
# Clean  #
##########
.PHONY: clean

clean: clean-executables
	rm -rf build/
	rm -f src/*.o
