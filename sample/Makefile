CC = g++

CXXFLAGS = -I../trimesh2/include
LDFLAGS = -g -L../trimesh2/lib.Linux64 -fopenmp
LDLIBS = -ltrimesh

all: model-info

model-info: model-info.o

clean:
	rm *.o
