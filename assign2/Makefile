CC = g++
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE = interp.x
SOURCES = interp.cpp linear.cpp spline.cpp triag_solve.cpp
OBJECTS = interp.o linear.o spline.o triag_solve.o

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(OBJECTS)

