CC=gcc
SRC=src
CFLAGS= -O3 -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wstrict-prototypes -Wwrite-strings -g

all: transverseHFK

transverseHFK: src/TransverseHFK.c
	$(CC) $(CFLAGS) $(SRC)/TransverseHFK.c -o transverseHFK

clean:
	rm -f transverseHFK

test: test_m10_132 test_m12n200 test_p_-4_-3_3 test_p1_-6_-3_3 test_p2_-6_-3_3 test_pr_-4_-3_3

test_long: test 

test_m10_132:
test_m12n200:
test_p_-4_-3_3:
test_p1_-6_-3_3:
test_p2_-6_-3_3:
test_pr_-4_-3_3:
