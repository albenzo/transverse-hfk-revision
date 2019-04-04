CC=gcc
SRC=src
CFLAGS= -O3 -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wstrict-prototypes -Wwrite-strings -g
LDFLAGS= 
TEST_DIR=test

UNAME_S= $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	CFLAGS += -I/usr/local/include/
	LDFLAGS += -L/usr/local/lib/ -largp
endif

ALL_TESTS := $(addsuffix .test, $(patsubst $(TEST_DIR)/%.in,%, $(wildcard $(TEST_DIR)/*.in)))

all: transverseHFK

python-install: src/TransverseHFK.c tHFK/__init__.py tHFK/tHFK.py tHFK/_transverseHFKmodule.c setup.py
	python setup.py install

install: transverseHFK
	cp ./transverseHFK /usr/bin/transverseHFK

uninstall:
	rm -rf /usr/bin/transverseHFK

transverseHFK: src/TransverseHFK.c src/states.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRC)/TransverseHFK.c -o transverseHFK

clean:
	rm -f transverseHFK
clean-python:
	python setup.py clean

test: $(ALL_TESTS)

%.test: transverseHFK $(TEST_DIR)/%.in $(TEST_DIR)/%.out
	./transverseHFK `cat $(TEST_DIR)/$*.in` 2>&1 | diff -q $(TEST_DIR)/$*.out - > /dev/null || (echo "Target $@ failed" && exit 1)

.PHONY: clean test %.test python-install clean-python install uninstall
