#
# OrderingSIMD
# by Miguel Cárcamo

INC_DIRS = -Isrc

CFLAGS = -c -w -Wall -msse4.1 

LDFLAGS = -fopenmp



main: build/simdsort.o 
	@ g++ $(LDFLAGS) build/*.o -o bin/simdsort.run
	@ echo "The compilation has been completed successfully"

build/simdsort.o: src/simdsort.cpp
	@ echo "src/simdsort.cpp"
	@ g++  $(CFLAGS) $(INC_DIRS) $(LDFLAGS) src/simdsort.cpp -o build/simdsort.o 

clean:
	@ clear
	@ echo "Cleaning folders.."
	@ rm -rf build/*
	@ rm -rf bin/*
	@ rm -rf output/*

run:
	@ clear
	@ echo "Execute Ordering SIMD"
	@ ./bin/simdsort.run -i ./input/100num.raw -o ./output/ordenada.raw -N 1024 -d 1 -L 4

simd: clean main run
