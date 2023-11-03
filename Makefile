CC = gcc
CFLAGS:= -std=c99 \
	-Wall \
	-Wextra \
	-pedantic \
	-Wfatal-errors

TEST_CFLAGS:=-std=c99

all: test
.PHONY: test
test:
	$(CC) $(TEST_CFLAGS) -I./modules/log/src -DLOG_USE_COLOR  -c ./modules/log/src/log.c -o .obj/log.o
	$(CC) $(TEST_CFLAGS) -I./modules/unity/src -DUNITY_INCLUDE_DOUBLE -c ./modules/unity/src/unity.c -o .obj/unity.o
	$(CC) $(TEST_CFLAGS) -I./src -c ./src/num/new.c -o .obj/new.o
	$(CC) $(TEST_CFLAGS) -I./src -I./modules/log/src -DM_PI=3.14159265359 -D_TOLERANCE=1e-10 -c ./src/num/num.c -o .obj/num.o
	$(CC) $(TEST_CFLAGS) -I./src/num -I./modules/log/src -DM_PI=3.14159265359 -c ./src/partition/partition.c -o .obj/partition.o
	$(CC) $(TEST_CFLAGS) -I./src/num -I./modules/log/src -DM_PI=3.14159265359 -c ./src/algorithm/integrate.c -o .obj/integrate.o
	$(CC) $(TEST_CFLAGS) -I./src/num -I./modules/log/src -DM_PI=3.14159265359 -c ./src/algorithm/algorithm.c -o .obj/algorithm.o
	$(CC) $(TEST_CFLAGS) -I./include -I./src/num -I./modules/log/src -DM_PI=3.14159265359 -c ./src/mittleff.c -o .obj/mittleff.o
	$(CC) $(TEST_CFLAGS) -I./include -I./modules/unity/src -I./modules/log/src  -DLOGLEVEL=1 -DUNITY_INCLUDE_DOUBLE -c ./test/test.c -o .obj/test.o
	$(CC) $$(ls .obj/*.o) -o test.out -lm -larb -lflint
	./test.out

.PHONY: clean clean-all
clean:
	@rm -rf .obj && mkdir .obj

