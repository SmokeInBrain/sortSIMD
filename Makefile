#
# OrderingSIMD
# by Miguel Cárcamo

INC_DIRS = -Isrc

CFLAGS = -c -w -Wall -msse4.1 

LDFLAGS = -fopenmp



main: build/simdsort.o build/minheap.o
	@ g++ $(LDFLAGS) build/*.o -o bin/simdsort.run
	@ echo "The compilation has been completed successfully"

build/simdsort.o: src/simdsort.cpp
	@ echo "src/simdsort.cpp"
	@ g++  $(CFLAGS) $(INC_DIRS) $(LDFLAGS) src/simdsort.cpp -o build/simdsort.o 

build/minheap.o: src/minheap.cpp
	@ echo "src/minheap.cpp"
	@ g++  $(CFLAGS) $(INC_DIRS) $(LDFLAGS) src/minheap.cpp -o build/minheap.o 

clean:
	@ clear
	@ echo "Cleaning folders.."
	@ rm -rf build/*
	@ rm -rf bin/*
	@ rm -rf output/*

run:
	@ clear
	@ echo "Execute Ordering SIMD"
	@ ./bin/simdsort.run -i ./input/16Mnum.raw -o ./output/ordenada.raw -N 16777216 -d 0 -L 2

simd: clean main run
