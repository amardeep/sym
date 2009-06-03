CC = g++

CXXFLAGS = -I../trimesh2/include -I.. 
LDFLAGS = -g -L../trimesh2/lib.Linux64 -fopenmp -L/usr/lib
LDLIBS = -ltrimesh -lGL -lGLU -lglut

all: sym-detect mesh_sym sym-view-curv mesh_sym1 plot-samples find-corr view-trans

sym-detect: sym-detect.o

clean:
	rm *.o
