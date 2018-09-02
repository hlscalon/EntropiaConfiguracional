main:
	g++ -O3 -Wall -shared -std=c++14 -Ilibs/Aboria/src -Ilibs/Aboria/third-party -Ilibs/pybind11/include -fPIC `python-config --includes` -DABORIA_LOG_LEVEL=1 *.cpp -o boost_graph.so
