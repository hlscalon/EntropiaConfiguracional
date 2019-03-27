CC = g++
CFLAGS = -O3 -Wall -shared -std=c++14
LIBS = -Isrc/cpp/libs/Aboria/src -Isrc/cpp/libs/Aboria/third-party -Isrc/cpp/libs/pybind11/include -Isrc/cpp/libs/nauty
OBJECTS = src/cpp/libs/nauty/nauty.a
PYTHON = -fPIC `python-config --includes`
OPENMP = -fopenmp
DEFINES = -DABORIA_LOG_LEVEL=1 -LOG #-DDEBUG
OBJECT = boost_graph.so
SOURCES = src/cpp/*.cpp
DEBUG_INFO = -g

main:
	$(CC) $(CFLAGS) $(LIBS) $(DEBUG_INFO) $(PYTHON) $(OPENMP) $(DEFINES) $(SOURCES) $(OBJECTS) -o $(OBJECT)
	mv $(OBJECT) src/python
