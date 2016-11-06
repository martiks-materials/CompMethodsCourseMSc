CC = g++
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE = machine.x
SOURCES = machine.cpp
OBJECTS = machine.o

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(OBJECTS)

