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
	$(CC) -I./modules/Unity/src -DUNITY_INCLUDE_DOUBLE -c ./modules/Unity/src/unity.c -o build_aux/unity.o
	$(CC) -I./src/num -c ./src/num/num.c -o build_aux/num.o
	$(CC) -I./modules/Unity/src -I./src/num -c test/test_num.c -o build_aux/test_num.o
	$(CC) build_aux/unity.o build_aux/num.o build_aux/test_num.o -o build/test_num.x -larb -lflint
	$(VALGRIND) --log-file=valgrind-test_num.txt ./build/test_num.x

# DEV_CFLAGS:=-std=c99
# DEV_INCDIR:=-I./include -I./src/include

# LOGLEVEL ?= 1
# all: test
# .PHONY: test
# test:
# 	@mkdir --parents build_aux
# 	$(CC) $(DEV_CFLAGS) -I./modules/log.c/src -DLOG_USE_COLOR -c ./modules/log.c/src/log.c -o build_aux/log.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/Unity/src -DUNITY_INCLUDE_DOUBLE -c ./modules/Unity/src/unity.c -o build_aux/unity.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -c ./modules/num.c/src/new.c -o build_aux/new.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -c ./modules/num.c/src/num.c -o build_aux/num.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/integration.c/include -c ./modules/integration.c/src/qsimp.c -o build_aux/qsimp.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/integration.c/include -c ./modules/integration.c/src/qtanhsinh.c -o build_aux/qtanhsinh.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/log.c/src -I./modules/num.c/include -I./modules/integration.c/include -I./modules/integration.c/src/include -c ./modules/integration.c/src/integration.c -o build_aux/integration.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/log.c/src -c $(DEV_INCDIR) -I./src -DM_PI=3.14159265359 -c ./src/partition.c -o build_aux/partition.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/integration.c/include -I./modules/log.c/src -c $(DEV_INCDIR) -I./src -DM_PI=3.14159265359 -c ./src/integrate.c -o build_aux/integrate.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/log.c/src -c $(DEV_INCDIR) -I./src -DM_PI=3.14159265359 -c ./src/algorithm.c -o build_aux/algorithm.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/log.c/src -c $(DEV_INCDIR) -I./src -DM_PI=3.14159265359 -c ./src/mittleff.c -o build_aux/mittleff.o
# 	$(CC) $(DEV_CFLAGS) $(DEV_INCDIR) \
# 	-I./modules/log.c/src \
# 	-I./modules/Unity/src  \
# 	-I./modules/num.c/include -I./module/log.c/src \
# 	-I./src -I./test  \
# 	$(DEV_INCDIR) \
# 	-DLOGLEVEL=$(LOGLEVEL) -DUNITY_INCLUDE_DOUBLE -c ./test/test.c -o build_aux/test.o
# 	$(CC) $$(ls build_aux/*.o) -o test.out -lm -larb -lflint -lgsl -lgslcblas
# 	@valgrind \
# 	--leak-check=full \
# 	--show-leak-kinds=all \
#         --track-origins=yes \
#         --verbose          \
#         --log-file=valgrind-out.txt \
#         ./test.out



# .PHONY: clean clean-all
# clean:
# 	@rm -rf build_aux && mkdir build_aux

