CC = g++ -std=c++11
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE1 = rosen.x 
EXECUTABLE2 = compare.x
EXECUTABLE3 = multivariate.x
SOURCES1 = rosen.cpp mcmc.cpp 
SOURCES2 = compare.cpp mcmc.cpp
SOURCES3 = multivariate.cpp mcmc.cpp
OBJECTS1 = $(SOURCES1:.cpp=.o)
OBJECTS2 = $(SOURCES2:.cpp=.o)
OBJECTS3 = $(SOURCES3:.cpp=.o)

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
	rm $(EXECUTABLE1) $(EXECUTABLE2) $(EXECUTABLE3) rosen.o compare.o multivariate.o mcmc.o

