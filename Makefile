# Makefile
#
TARGET_J  = bin/poisson_j        # Jacobi
TARGET_GS = bin/poisson_gs       # Gauss-Seidel

SOURCES = main.c print.c alloc3d.c define_u_f.c
OBJECTS = bin/print.o bin/alloc3d.o bin/define_u_f.o
MAIN_J  = bin/main_j.o
MAIN_GS = bin/main_gs.o
OBJS_J  = $(MAIN_J) bin/jacobi.o
OBJS_GS = $(MAIN_GS) bin/gauss_seidel.o
TEST_SRC = test.c
TEST_OUT = test

# options and settings for the GCC compilers
#
CC      = gcc
DEFS    =
OPT     = -g
IPO     =
ISA     =
CHIP    =
ARCH    =
PARA    =
CFLAGS  = $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) $(XOPTS)
LDFLAGS = -lm

# Ensure bin directory exists before compiling
all: bin $(TARGET_J) $(TARGET_GS)

bin:
	@mkdir -p bin

$(TARGET_J): $(OBJECTS) $(OBJS_J)
	$(CC) -o $@ $(CFLAGS) $(OBJS_J) $(OBJECTS) $(LDFLAGS)

$(TARGET_GS): $(OBJECTS) $(OBJS_GS)
	$(CC) -o $@ $(CFLAGS) $(OBJS_GS) $(OBJECTS) $(LDFLAGS)

bin/main_j.o: main.c
	$(CC) -o $@ -D_JACOBI $(CFLAGS) -c main.c

bin/main_gs.o: main.c
	$(CC) -o $@ -D_GAUSS_SEIDEL $(CFLAGS) -c main.c

bin/%.o: %.c
	$(CC) -o $@ $(CFLAGS) -c $<

test:
	$(CC) $(CFLAGS) $(TEST_SRC) $(LDFLAGS) && ./a.out
	@rm -f a.out

clean:
	@/bin/rm -f core bin/*.o *~

realclean: clean
	@/bin/rm -f bin/poisson_j bin/poisson_gs
	@rmdir bin || true

# DO NOT DELETE

bin/main_j.o: main.c print.h jacobi.h
bin/main_gs.o: main.c print.h gauss_seidel.h
bin/print.o: print.h
bin/define_u_f.o: define_u_f.c define_u_f.h