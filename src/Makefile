CC = g++

CXXFLAGS = -I../trimesh2/include -I.. 
LDFLAGS = -g -L../trimesh2/lib.Linux64 -fopenmp -L/usr/lib
LDLIBS = -ltrimesh -lGL -lGLU -lglut

all: sym-detect

sym-detect: sym-detect.o

clean:
	rm *.o
