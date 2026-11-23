.PHONY: clean test debug all
CC=gcc
CXX=g++

all:
	cmake -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
	cmake --build build
	mv build/dist .
	mv build/dec .
	mv build/part .
	mv build/toperm .

debug:
	cmake -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX)  -H. -Bbuild -DCMAKE_BUILD_TYPE=Debug
	cmake --build build
	mv build/dist .
	mv build/dec .
	mv build/part .
	mv build/toperm .

test:
	cmake -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX)  -H. -Bbuild -DCMAKE_BUILD_TYPE=Debug
	cmake --build build
	cd build && ctest --output-on-failure

clean:
	cd build && make clean
