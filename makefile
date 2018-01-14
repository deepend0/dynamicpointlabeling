all:
	g++ -g -O0 -m64 -std=c++14 -o bin/test -m64 src/*.cpp src/test/*.cpp -Isrc -Isrc/test -I"/home/oakile/Workspace/lib/galib247" -L"/home/oakile/Workspace/lib/galib247/ga" -lga -I"/home/oakile/Workspace/lib/boost_1_60_0/" -L"/home/oakile/Workspace/lib/boost_1_60_0/stage/lib" -lboost_filesystem -lboost_system -lboost_regex
	
clean:	
	rm bin/test