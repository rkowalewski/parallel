CC=gcc-4.9
FLAGS=-std=c99
MPICC=mpicc

all :
	gcc ${FLAGS} -o filter filter.c
	${MPICC} ${FLAGS} -o filter_p filter_p.c

clean : 
	rm -f *~
	rm -f a.out
	rm -f *.ppm
	rm -f *.pgm
	rm -f filter
	rm -f filter_p

