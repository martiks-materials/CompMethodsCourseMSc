CC = g++
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE1 = uniform.x 
EXECUTABLE2 =transform.x 
EXECUTABLE3 = reject.x
SOURCES1 = uniform.cpp 
SOURCES2 = transform.cpp 
SOURCES3 = reject.cpp compare.cpp
OBJECTS1 = uniform.o 
OBJECTS2 = transform.o 
OBJECTS3 = reject.o compare.o

all: $(SOURCES1) $(SOURCES2) $(SOURCES3) $(EXECUTABLE1) $(EXECUTABLE2) $(EXECUTABLE3)


$(EXECUTABLE1): $(OBJECTS1)
	$(CC) $(OBJECTS1) -o $@ $(LDFLAGS)

$(EXECUTABLE2): $(OBJECTS2)
	$(CC) $(OBJECTS2) -o $@ $(LDFLAGS)

$(EXECUTABLE3): $(OBJECTS3)
	$(CC) $(OBJECTS3) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE1) $(EXECUTABLE2) $(EXECUTABLE3) $(OBJECTS1) $(OBJECTS2) $(OBJECTS3)

