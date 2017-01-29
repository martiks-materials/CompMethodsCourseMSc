CC = g++ -std=c++11
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE1 = erfrand.x 
SOURCES1 = erfrand.cpp mcint.cpp multifunc.cpp
OBJECTS1 = $(SOURCES1:.cpp=.o)

all: $(SOURCES1) $(EXECUTABLE1)

$(EXECUTABLE1): $(OBJECTS1)
	$(CC) $(OBJECTS1) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE1) $(OBJECTS1)

