CC=gcc
SRC=src
CFLAGS= -O3 -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wstrict-prototypes -Wwrite-strings -g

all: transverseHFK

transverseHFK: src/TransverseHFK.c
	$(CC) $(CFLAGS) $(SRC)/TransverseHFK.c -o transverseHFK

clean:
	rm -f transverseHFK
