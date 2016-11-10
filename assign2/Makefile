CC = g++
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE = interp.x
SOURCES = interp.cpp linear.cpp triag_own.cpp spline.cpp
OBJECTS = interp.o linear.o triag_own.o spline.o

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(OBJECTS)

