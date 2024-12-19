a:integrate.o main.o
	g++ -std=c++11 integrate.o main.o
	
main.o:main.cpp integrate.h
	g++ -std=c++11  main.cpp -c

integrate.o:integrate.cpp integrate.h
	g++ -std=c++11  integrate.cpp -c
