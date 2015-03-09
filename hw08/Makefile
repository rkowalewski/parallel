MPICC=mpicc
CFLAGS=-std=c99 -Wall -c
LDFLAGS=-Wall
SOURCES=stencil_testbench.c stencil.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTEABLE=stencil

all: $(EXECUTEABLE)

$(EXECUTEABLE): $(OBJECTS)
	$(MPICC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:	$(SOURCES)
	$(MPICC) $(CFLAGS) $< -o $@

clean : 
	rm -f *.pgm *.o $(EXECUTEABLE)

