CC = g++
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE = reject.x
SOURCES = reject.cpp compare.cpp
OBJECTS = reject.o compare.o

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(OBJECTS)

