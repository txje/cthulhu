CC=gcc
CFLAGS=-std=c99 -O2

OBJECTS = cthulhu covstat

all: $(OBJECTS)

cthulhu: src/main.c src/taxonomy.c
	$(CC) $(CFLAGS) src/main.c src/taxonomy.c src/paf.c incl/minimap2/libminimap2.a -o cthulhu -lhts -lz -lm -lpthread

covstat: src/covstat.c
	$(CC) $(CFLAGS) src/covstat.c -o covstat -lhts -lz -lm

.PHONY: clean
clean:
	-rm $(OBJECTS)
