CC=g++ -std=c++11 -Wall
PYTHON_LIB=`python3-config --includes` 
PYTHON_EXT=`python3-config --extension-suffix`

default: All

All: similarity.cpp
	$(CC) -fopenmp -O3 -shared -fPIC $(PYTHON_LIB) -I `python3 -m pybind11 --includes` -c similarity.cpp
	$(CC) -fopenmp -O3 -shared -fPIC $(PYTHON_LIB) -I `python3 -m pybind11 --includes` similarity.o -o similarity$(PYTHON_EXT)

clean:
	rm -f *.o
