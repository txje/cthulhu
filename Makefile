CC=gcc
CFLAGS=-std=c99 -O2

OBJECTS = cattax covstat

all: $(OBJECTS)

cattax: src/cattax.c src/taxonomy.c
	$(CC) $(CFLAGS) src/cattax.c src/taxonomy.c incl/minimap2/libminimap2.a -o cattax -lhts -lz -lm -lpthread

covstat: src/covstat.c
	$(CC) $(CFLAGS) src/covstat.c -o covstat -lhts -lz -lm

.PHONY: clean
clean:
	-rm $(OBJECTS)
