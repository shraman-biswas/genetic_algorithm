CXX	=	gcc
FLAGS	=	-Wall -O3 -ggdb
SOURCE	=	main.c genetic_algorithm.c
BIN	=	main

build:
	$(CXX) $(FLAGS) $(SOURCE) -o $(BIN)

clean:
	rm -f *~ $(BIN)
