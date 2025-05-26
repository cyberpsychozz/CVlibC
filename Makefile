# Makefile for CVlibC project

CC = gcc
CFLAGS = -Wall -O2
LDFLAGS = -lm

SRC = user_friendly_main.c
TARGET = cvlib_demo

# Default rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Clean rule
clean:
	rm -f $(TARGET)
