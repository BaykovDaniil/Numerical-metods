a:derivative.o main.o
	g++ -std=c++11 derivative.o main.o
	
main.o:main.cpp derivative.h
	g++ -std=c++11  main.cpp -c

derivative.o:derivative.cpp derivative.h
	g++ -std=c++11  derivative.cpp -c
