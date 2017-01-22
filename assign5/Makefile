CC = g++ -std=c++11
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE1 = intmain.x
EXECUTABLE2 = erfmain.x 
SOURCES1 = intmain.cpp rk4.cpp specifics.cpp 
SOURCES2 = erfmain.cpp rk4.cpp
OBJECTS1 = $(SOURCES1:.cpp=.o) 
OBJECTS2 = $(SOURCES2:.cpp=.o)

all: $(SOURCES1) $(SOURCES2) $(EXECUTABLE1) $(EXECUTABLE2)

$(EXECUTABLE1): $(OBJECTS1)
	$(CC) $(OBJECTS1) -o $@ $(LDFLAGS)

$(EXECUTABLE2): $(OBJECTS2)
	$(CC) $(OBJECTS2) -o $@ $(LDFLAGS)
	
%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE1) $(OBJECTS1) $(EXECUTABLE2) $(OBJECTS2)

