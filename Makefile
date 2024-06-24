CC = gcc

CFLAGS = -Wall -Wextra -std=c11 -Iinclude

TARGET = build/w3
SRCDIR = src
INCLUDEDIR = include

SRCS = $(SRCDIR)/w3.c $(SRCDIR)/input.c $(SRCDIR)/equilibrium.c
OBJS = $(SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

src/w3.o: src/w3.c include/input.h include/equilibrium.h
	$(CC) $(CFLAGS) -c $< -o $@

src/input.o: src/input.c include/input.h
	$(CC) $(CFLAGS) -c $< -o $@

src/equilibrium.o: src/equilibrium.c include/equilibrium.h
	$(CC) $(CFLAGS) -c $< -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(SRCDIR)/*.o $(TARGET)

.PHONY: all clean		
