APP_NAME=graph
APP_NAME2=validateMST

SRCS = graph.cpp
OBJS = $(SRCS:.cpp=.o)

SRCS2 = validateMST.cpp
OBJS2 = $(SRCS2:.cpp=.o)

CXX = g++
CXXFLAGS = -Wall -O3 -std=c++17 -m64 -I. -fopenmp -Wno-unknown-pragmas -Wno-maybe-uninitialized

all: $(APP_NAME) $(APP_NAME2)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

$(APP_NAME2): $(OBJS2)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS2)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) $(APP_NAME2) *.class