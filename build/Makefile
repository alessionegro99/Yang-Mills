CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99 

binaries = u1_conf

all: $(binaries)

boxmuller.o: ../lib/boxmuller.c ../include/boxmuller.h
	$(CC) $(CFLAGS) -c ../lib/boxmuller.c

geometry.o: ../lib/geometry.c ../include/geometry.h 
	$(CC) $(CFLAGS) -c ../lib/geometry.c  

measure.o: ../lib/measure.c ../include/measure.h
	$(CC) $(CFLAGS) -c ../lib/measure.c

pcg32min.o: ../lib/pcg32min.c ../include/pcg32min.h 
	$(CC) $(CFLAGS) -c ../lib/pcg32min.c

random.o: ../lib/random.c ../include/pcg32min.h ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

read_data.o: ../lib/read_data.c ../include/read_data.h
	$(CC) $(CFLAGS) -c ../lib/read_data.c

update.o: ../lib/update.c ../include/update.h
	$(CC) $(CFLAGS) -c ../lib/update.c

u1_conf: ../src/u1_conf.c geometry.o random.o pcg32min.o measure.o update.o
	$(CC) $(CFLAGS) -c ../src/u1_conf.c  
	$(CC) $(CFLAGS) u1_conf.o geometry.o random.o pcg32min.o measure.o update.o -o $@ -lm

.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o
