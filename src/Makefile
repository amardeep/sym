CC = g++

CXXFLAGS = -I../trimesh2/include -I.. -I../ann/include/ANN
LDFLAGS = -g -L../trimesh2/lib.Linux64 -fopenmp -L/usr/lib -L../ann/lib
LDLIBS = -ltrimesh -lGL -lGLU -lglut -lANN

all: sym-detect mesh_sym sym-view-curv mesh_sym1 plot-samples find-corr view-trans sym-view-curv-kd mesh_sym2

sym-detect: sym-detect.o

clean:
	rm *.o
