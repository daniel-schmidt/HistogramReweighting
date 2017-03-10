TARGET = multihist
LIBS = -lgsl -lgslcblas -lm
CC = clang-3.6
CFLAGS = -std=gnu11 -Wall -g -O3 -march=native

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@
clean:
	-rm -f *.o
	-rm -f $(TARGET)