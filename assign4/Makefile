CC = g++ -std=c++11
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE1 = numint.x 
SOURCES1 = numint.cpp  
OBJECTS1 = numint.o

all: $(SOURCES1) $(EXECUTABLE1)

$(EXECUTABLE1): $(OBJECTS1)
	$(CC) $(OBJECTS1) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE1) $(OBJECTS1)

