CC = g++
CFLAGS = -O3 -Wall -shared -std=c++14
LIBS = -Isrc/cpp/libs/Aboria/src -Isrc/cpp/libs/Aboria/third-party -Isrc/cpp/libs/pybind11/include -Isrc/cpp/libs/nauty
OBJECTS = src/cpp/libs/nauty/nauty.a
PYTHON = -fPIC `python-config --includes`
OPENMP = -fopenmp
DEFINES = -DABORIA_LOG_LEVEL=1 -LOG #-DDEBUG
SOURCES = src/cpp/*.cpp
# SOURCES_OUT = src/cpp/FindIsomorphicIndex.cpp
# SOURCES = $(filter-out $(SOURCES_OUT), $(wildcard src/cpp/*.cpp))
# SOURCES_GPU = src/cpp/IsIsomorphic.cu src/cpp/FindIsomorphicIndex.cu

OBJECT = boost_graph.so

main:
	# $(CC) $(CFLAGS) $(LIBS) $(PYTHON) $(OPENMP) $(DEFINES) $(SOURCES) $(SOURCES_OUT) $(OBJECTS) -o $(OBJECT)
	$(CC) $(CFLAGS) $(LIBS) $(PYTHON) $(OPENMP) $(DEFINES) $(SOURCES) $(OBJECTS) -o $(OBJECT)
	mv $(OBJECT) src/python

gccold:
	scl run devtoolset-7 bash
	gcc -v

gpu:
	nvcc -c src/cpp/FindIsomorphicIndex.cu -o findIso.o -std=c++14 --expt-extended-lambda
	$(CC) $(CFLAGS) $(LIBS) $(PYTHON) $(DEFINES) $(SOURCES) $(OBJECTS) findIso.o -o $(OBJECT)
	mv $(OBJECT) src/python
