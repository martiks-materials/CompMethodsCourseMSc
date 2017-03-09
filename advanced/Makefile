CC = g++ -std=c++11
CFLAGS = -c -Wall -Wextra
LDFLAGS =

EXECUTABLE1 = motion.x
EXECUTABLE2 = fourtet.x
EXECUTABLE3 = doubslit.x
SOURCES1 = motion.cpp fft_gen.cpp
SOURCES2 = fourtet.cpp fft_gen.cpp func_zoo.cpp
SOURCES3 = doubslit.cpp fft_gen.cpp
OBJECTS1 = $(SOURCES1:.cpp=.o) 
OBJECTS2 = $(SOURCES2:.cpp=.o) 
OBJECTS3 = $(SOURCES3:.cpp=.o) 

all: $(SOURCES1) $(EXECUTABLE1) $(SOURCES2) $(EXECUTABLE2) $(SOURCES3) $(EXECUTABLE3)

$(EXECUTABLE1): $(OBJECTS1)
	$(CC) $(OBJECTS1) -o $@ $(LDFLAGS)

$(EXECUTABLE2): $(OBJECTS2)
	$(CC) $(OBJECTS2) -o $@ $(LDFLAGS)

$(EXECUTABLE3): $(OBJECTS3)
	$(CC) $(OBJECTS3) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm  $(EXECUTABLE1) $(OBJECTS1) $(EXECUTABLE2) $(OBJECTS2) $(EXECUTABLE3) $(OBJECTS3)

