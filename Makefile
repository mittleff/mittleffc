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
	$(CC) $(CFLAGS) -I./modules/Unity/src -DUNITY_INCLUDE_DOUBLE -c ./modules/Unity/src/unity.c -o build/unity.o
	$(CC) $(CFLAGS) -I./src/num -c ./src/num/new.c -o build/new.o
	$(CC) $(CFLAGS) -I./src/num -c ./src/num/num.c -o build/num.o
	$(CC) $(CFLAGS) -I./modules/Unity/src -I./src/num -DUNITY_INCLUDE_DOUBLE -c test/test_num.c -o build/test_num.o
	$(CC) $(CFLAGS) build/unity.o build/new.o build/num.o build/test_num.o -o build/test_num.x -lflint -lgsl -lgslcblas -lm
	$(VALGRIND) --log-file=valgrind-test_num.txt ./build/test_num.x

test_mittleff:
	@mkdir --parents build
	$(CC) $(CFLAGS) -I./modules/Unity/src -DUNITY_INCLUDE_DOUBLE -c ./modules/Unity/src/unity.c -o build/unity.o
	$(CC) $(CFLAGS) -I./modules/log.c/src -DLOG_USE_COLOR -c ./modules/log.c/src/log.c -o build/log.o
	$(CC) $(CFLAGS) -I./modules/log.c/src -I./src/num -c ./src/num/num.c -o build/num.o
	$(CC) $(CFLAGS) -I./modules/log.c/src -I./include -I./src/include -I./src/num -c ./src/partition.c -o build/partition.o
	$(CC) $(CFLAGS) -I./modules/log.c/src -I./include -I./src/include -I./src/num -c ./src/algorithm.c -o build/algorithm.o
	$(CC) $(CFLAGS) -I./modules/log.c/src -I./include -I./src/include -I./src/num -c ./src/quad.c -o build/quad.o
	$(CC) $(CFLAGS) -I./modules/log.c/src -I./include -I./src/include -I./src/num -c ./src/integrate.c -o build/integrate.o
	$(CC) $(CFLAGS) -I./modules/log.c/src -I./include -I./src/include -I./src/num -c ./src/mittleff.c -o build/mittleff.o
	$(CC) $(CFLAGS) -I./include -I./modules/Unity/src -I./src/num -DUNITY_INCLUDE_DOUBLE -c test/test_mittleff.c -o build/test_mittleff.o
	$(CC) $(CFLAGS) build/unity.o build/log.o build/algorithm.o build/quad.o build/integrate.o build/num.o build/partition.o build/mittleff.o build/test_mittleff.o -o build/test_mittleff.x -lflint -lgsl -lgslcblas -lm
	$(VALGRIND) --log-file=valgrind-test_mittleff.txt ./build/test_mittleff.x

# DEV_CFLAGS:=-std=c99
# DEV_INCDIR:=-I./include -I./src/include

# LOGLEVEL ?= 1
# all: test
# .PHONY: test
# test:
# 	@mkdir --parents build
# 	$(CC) $(DEV_CFLAGS) -I./modules/log.c/src -DLOG_USE_COLOR -c ./modules/log.c/src/log.c -o build/log.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/Unity/src -DUNITY_INCLUDE_DOUBLE -c ./modules/Unity/src/unity.c -o build/unity.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -c ./modules/num.c/src/new.c -o build/new.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -c ./modules/num.c/src/num.c -o build/num.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/integration.c/include -c ./modules/integration.c/src/qsimp.c -o build/qsimp.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/integration.c/include -c ./modules/integration.c/src/qtanhsinh.c -o build/qtanhsinh.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/log.c/src -I./modules/num.c/include -I./modules/integration.c/include -I./modules/integration.c/src/include -c ./modules/integration.c/src/integration.c -o build/integration.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/log.c/src -c $(DEV_INCDIR) -I./src -DM_PI=3.14159265359 -c ./src/partition.c -o build/partition.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/integration.c/include -I./modules/log.c/src -c $(DEV_INCDIR) -I./src -DM_PI=3.14159265359 -c ./src/integrate.c -o build/integrate.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/log.c/src -c $(DEV_INCDIR) -I./src -DM_PI=3.14159265359 -c ./src/algorithm.c -o build/algorithm.o
# 	$(CC) $(DEV_CFLAGS) -I./modules/num.c/include -I./modules/log.c/src -c $(DEV_INCDIR) -I./src -DM_PI=3.14159265359 -c ./src/mittleff.c -o build/mittleff.o
# 	$(CC) $(DEV_CFLAGS) $(DEV_INCDIR) \
# 	-I./modules/log.c/src \
# 	-I./modules/Unity/src  \
# 	-I./modules/num.c/include -I./module/log.c/src \
# 	-I./src -I./test  \
# 	$(DEV_INCDIR) \
# 	-DLOGLEVEL=$(LOGLEVEL) -DUNITY_INCLUDE_DOUBLE -c ./test/test.c -o build/test.o
# 	$(CC) $$(ls build/*.o) -o test.out -lm -larb -lflint -lgsl -lgslcblas
# 	@valgrind \
# 	--leak-check=full \
# 	--show-leak-kinds=all \
#         --track-origins=yes \
#         --verbose          \
#         --log-file=valgrind-out.txt \
#         ./test.out



# .PHONY: clean clean-all
# clean:
# 	@rm -rf build && mkdir build

