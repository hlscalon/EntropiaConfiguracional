CC = g++
CFLAGS = -O3 -Wall -shared -std=c++14
LIBS = -Isrc/cpp/libs/Aboria/src -Isrc/cpp/libs/Aboria/third-party -Isrc/cpp/libs/pybind11/include
PYTHON = -fPIC `python-config --includes`
OPENMP = -fopenmp
DEFINES = -DABORIA_LOG_LEVEL=1
SOURCES_OUT = src/cpp/FindIsomorphicIndex.cpp src/cpp/IsIsomorphic.cpp
SOURCES = $(filter-out $(SOURCES_OUT), $(wildcard src/cpp/*.cpp))
#SOURCES_GPU = src/cpp/*.cu
SOURCES_GPU = src/cpp/IsIsomorphic.cu src/cpp/FindIsomorphicIndex.cu

OBJECT = boost_graph.so

main:
	$(CC) $(CFLAGS) $(LIBS) $(PYTHON) $(OPENMP) $(DEFINES) $(SOURCES) $(SOURCES_OUT) -o $(OBJECT)
	mv $(OBJECT) src/python

gccold:
	scl run devtoolset-7 bash
	gcc -v

gpu:
	nvcc -c src/cpp/IsIsomorphic.cu -o isIso.o -std=c++14 --expt-extended-lambda
	nvcc -c src/cpp/FindIsomorphicIndex.cu -o findIso.o -std=c++14 --expt-extended-lambda
	$(CC) $(CFLAGS) $(LIBS) $(PYTHON) $(DEFINES) $(SOURCES) isIso.o findIso.o -o $(OBJECT)
	mv $(OBJECT) src/python