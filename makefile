CC=gcc
SRC=src
CFLAGS= -O3 -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wstrict-prototypes -Wwrite-strings -g

all: transverseHFK

transverseHFK: src/TransverseHFK.c
	$(CC) $(CFLAGS) $(SRC)/TransverseHFK.c -o transverseHFK

clean:
	rm -f transverseHFK

test: test_m10_132 test_10_132 test_1_m12n200 test_2_m12n200 test_p1_-4_-3_3 test_p2_-4_-3_3 test_p1_-6_-3_3 test_p2_-6_-3_3

test_long: test 

test_m10_132:
	./transverseHFK `cat test/m10_132.in)` | diff -q test/m10_132.out - > /dev/null || (echo "Target $@ failed" && exit 1)
test_10_132:
	./transverseHFK `cat test/10_132.in)` | diff -q test/10_132.out - > /dev/null || (echo "Target $@ failed" && exit 1)
test_1_m12n200:
	./transverseHFK `cat test/1_m12n200.in)` | diff -q test/1_m12n200.out - > /dev/null || (echo "Target $@ failed" && exit 1)
test_2_m12n200:
	./transverseHFK `cat test/2_m12n200.in)` | diff -q test/2_m12n200.out - > /dev/null || (echo "Target $@ failed" && exit 1)
test_p1_-4_-3_3:
	./transverseHFK `cat test/p1_-4_-3_3.in)` | diff -q test/p1_-4_-3_3.out - > /dev/null || (echo "Target $@ failed" && exit 1)
test_p2_-4_-3_3:
	./transverseHFK `cat test/p2_-4_-3_3.in)` | diff -q test/p2_-4_-3_3.out - > /dev/null || (echo "Target $@ failed" && exit 1)
test_p1_-6_-3_3:
	./transverseHFK `cat test/p1_-6_-3_3.in` | diff -q test/p1_-6_-3_3.out - > /dev/null || (echo "Target $@ failed" && exit 1)
test_p2_-6_-3_3:
	./transverseHFK `cat test/p2_-6_-3_3.in` | diff -q test/p2_-6_-3_3.out - > /dev/null || (echo "Target $@ failed" && exit 1)
