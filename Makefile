CC = g++
CFLAGS = -O3 -Wall -shared -std=c++14
LIBS = -Isrc/cpp/libs/Aboria/src -Isrc/cpp/libs/Aboria/third-party -Isrc/cpp/libs/pybind11/include
PYTHON = -fPIC `python-config --includes`
OPENMP = -fopenmp
DEFINES = -DABORIA_LOG_LEVEL=1 #-DDEBUG
SOURCES = src/cpp/*.cpp
OBJECT = boost_graph.so

main:
	$(CC) $(CFLAGS) $(LIBS) $(PYTHON) $(OPENMP) $(DEFINES) $(SOURCES) -o $(OBJECT)
	mv $(OBJECT) src/python
