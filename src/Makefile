# Profiler
# PROFILER=-pg
PROFILER=

CFLAGS=-mavx2 -Wall -O2 -c $(PROFILER)

LDFLAGS=-lm $(PROFILER)

all: main

main: polynome.o main.o vect.o
	gcc polynome.o vect.o main.o $(LDFLAGS) -o main

polynome.o: polynome.h polynome.c
	gcc $(CFLAGS) polynome.c

vect.o: vect.h vect.c
	gcc $(CFLAGS) vect.c

main.o: main.c polynome.h
	gcc $(CFLAGS) main.c

clean: 
	rm -f *.o
	rm -f *.ps
	rm -f main
