CC = gcc
CFLAGS = -I./include -Wall -Wextra -std=c11

TARGET = build/w3


SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c, build/%.o, $(SRCS))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

build/%.o: src/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf build

.PHONY: all clean
