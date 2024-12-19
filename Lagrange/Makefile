a:lagrange.o main.o
	g++ -std=c++11 lagrange.o main.o
	
main.o:main.cpp lagrange.h
	g++ -std=c++11  main.cpp -c

lagrange.o:lagrange.cpp lagrange.h
	g++ -std=c++11  lagrange.cpp -c
