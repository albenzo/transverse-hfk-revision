CC=gcc
CFLAGS= -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wstrict-prototypes -Wwrite-strings -g -O3
LDFLAGS=
LIBS=
INCLUDES=
SRC_DIR=./src
PY_DIR=./tHFK
TEST_DIR=test
BUILD_DIR=./build
EXEC=transverseHFK

SRCS= $(shell find $(SRC_DIR) -name *.c)
OBJS= $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS= $(OBJS:.o=.d)

UNAME_S= $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	INCLUDES += -I/usr/local/include/
	LDFLAGS += -L/usr/local/lib/
	LIBS += -largp
endif

ALL_TESTS := $(addsuffix .test, $(patsubst $(TEST_DIR)/%.in,%, $(wildcard $(TEST_DIR)/*.in)))

all: $(BUILD_DIR)/$(EXEC)

python-install: $(SRC_DIR)/TransverseHFK.c $(PY_DIR)/__init__.py $(PY_DIR)/tHFK.py $(PY_DIR)/_transverseHFKmodule.c setup.py
	python setup.py install

install: $(EXEC)
	cp ./$(EXEC) /usr/bin/$(EXEC)

uninstall:
	rm /usr/bin/$(EXEC)

$(BUILD_DIR)/$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) $(LIBS) $(OBJS) -o $@

$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) $(LIBS) -c $< -o $@

clean: clean-python
	$(RM) -rf $(BUILD_DIR)

clean-python:
	python setup.py clean

test: $(ALL_TESTS)

%.test: $(BUILD_DIR)/$(EXEC) $(TEST_DIR)/%.in $(TEST_DIR)/%.out
	./$(BUILD_DIR)/$(EXEC) `cat $(TEST_DIR)/$*.in` 2>&1 | diff -q $(TEST_DIR)/$*.out - > /dev/null || (echo "Target $@ failed" && exit 1)

.PHONY: clean test %.test python-install clean-python install uninstall

-include $(DEPS)

MKDIR_P = mkdir -p
