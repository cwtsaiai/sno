CXX	 := g++ 
CXXFLAGS := -Wall -O3 -std=c++20
SRCS	 := main_snopso.cpp
OBJS	 := $(SRCS:.cpp=.o)
APPS	 := $(SRCS:.cpp=)
INCS	 := sno-de.h test_function.h

all: $(APPS)

$(APPS): $(OBJS)
	$(CXX) -o $@ $<

%.o: %.cpp $(INCS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	-rm -f $(APPS) *.o
