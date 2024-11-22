
CXX = g++
CFLAGS = -Wall -O2 -std=c++11 -msse4.2 -mbmi2
LDFLAGS = -msse4.2 -mbmi2


SRCS = test.cc crc.cpp

OBJS = $(SRCS:.cpp=.o)

EXE = test


all: $(EXE)


$(EXE): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS)


%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@


clean:
	rm -f $(EXE) $(OBJS)


.PHONY: all clean
