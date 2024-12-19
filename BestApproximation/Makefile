a:best_approx.o main.o
	g++ -std=c++11 best_approx.o main.o
	
main.o:main.cpp best_approx.h
	g++ -std=c++11  main.cpp -c

best_approx.o:best_approx.cpp best_approx.h
	g++ -std=c++11  best_approx.cpp -c
