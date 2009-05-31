CC = g++

CXXFLAGS = -I../trimesh2/include -I.. 
LDFLAGS = -g -L../trimesh2/lib.Linux64 -fopenmp 
LDLIBS = -ltrimesh

all: model-info curvature

model-info: model-info.o

clean:
	rm *.o
