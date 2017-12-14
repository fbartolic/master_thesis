CC 			= gcc 
MAKE 		= /usr/bin/make
CFLAGS 		= -c -Wall -DLINUX 
LIBDIRS 	= -L/usr/lib 
INCDIRS 	= -I/home/fran/Drive/Work/master_thesis/binary_c/src 
LIBS 		= -lbinary_c -lbfd -lm -lc 

all: evolve.o evolve.so clean

evolve.o: evolve.c
	$(CC)$ $(CFLAGS)$ evolve.c -o evolve.o $(LIBDIRS)$ $(LIBS)$ $(INCDIRS)$

evolve.so: evolve.o 
	$(CC)$	evolve.o -shared -o libevolve.so $(LIBDIRS)$ $(LIBS)$ 

clean:
	rm -f *.o 
