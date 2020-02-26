all:
	g++ -g -O0 -m64 -std=c++14 -o bin/test \
		src/main/graph/*.cpp \
		src/main/plp/*.cpp \
		src/main/gmxga/*.cpp \
		src/main/gaplp/*.cpp \
		src/main/simulation/*.cpp \
		src/main/simulation/confgraphgen/*.cpp \
		src/main/fhplp/*.cpp \
		src/test/*.cpp \
		-Isrc/main \
		-Isrc/main/graph \
		-Isrc/main/plp \
		-Isrc/main/gmxga \
		-Isrc/main/gaplp \
		-Isrc/main/simulation/ \
		-Isrc/main/simulation/confgraphgen \
		-Isrc/main/fhplp/ \
		-Isrc/test \
		-I"../FastStaticPointLabelPlacement/src" \
		-L"../FastStaticPointLabelPlacement/lib" -lfastheuristicplp \
		-I"../lib/galib247" \
		-L"../lib/galib247/ga" -lga \
		-I"../lib/boost_1_60_0/" \
		-L"../lib/boost_1_60_0/stage/lib" -lboost_filesystem -lboost_system -lboost_regex
	
clean:
	rm bin/test
