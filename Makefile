CC = gcc
CFLAGS:= -std=c99 \
	-Wall \
	-Wextra \
	-pedantic \
	-Wfatal-errors

DEV_CFLAGS:=-std=c99
DEV_INCDIR:=-I./include -I./modules/log/src 

all: test
.PHONY: test
test:
	$(CC) $(DEV_CFLAGS) -I./modules/log/src -DLOG_USE_COLOR  -c ./modules/log/src/log.c -o .obj/log.o
	$(CC) $(DEV_CFLAGS) -I./modules/unity/src -DUNITY_INCLUDE_DOUBLE -c ./modules/unity/src/unity.c -o .obj/unity.o
	$(CC) $(DEV_CFLAGS) $(DEV_INCDIR) -c ./src/new.c -o .obj/new.o
	$(CC) $(DEV_CFLAGS) $(DEV_INCDIR) -DM_PI=3.14159265359 -D_TOLERANCE=1e-10 -c ./src/num.c -o .obj/num.o
	$(CC) $(DEV_CFLAGS) $(DEV_INCDIR) -DM_PI=3.14159265359 -c ./src/partition.c -o .obj/partition.o
	$(CC) $(DEV_CFLAGS) $(DEV_INCDIR) -DM_PI=3.14159265359 -c ./src/integrate.c -o .obj/integrate.o
	$(CC) $(DEV_CFLAGS) $(DEV_INCDIR) -DM_PI=3.14159265359 -c ./src/algorithm.c -o .obj/algorithm.o
	$(CC) $(DEV_CFLAGS) $(DEV_INCDIR) -DM_PI=3.14159265359 -c ./src/mittleff.c -o .obj/mittleff.o
	$(CC) $(DEV_CFLAGS) $(DEV_INCDIR) -I./modules/unity/src -I./modules/log/src  -DLOGLEVEL=1 -DUNITY_INCLUDE_DOUBLE -c ./test/test.c -o .obj/test.o
	$(CC) $$(ls .obj/*.o) -o test.out -lm -larb -lflint
	./test.out

.PHONY: clean clean-all
clean:
	@rm -rf .obj && mkdir .obj

