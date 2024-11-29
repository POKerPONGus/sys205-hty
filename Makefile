
CPPFLAGS = -std=c++20 -g 

convert: src/csv_to_hty.cpp
	g++ $(CPPFLAGS) \
		-o bin/convert.out \
		src/csv_to_hty.cpp;

analyze: src/analyze.cpp
	g++ $(CPPFLAGS) \
		-o bin/analyze.out \
		src/analyze.cpp;

