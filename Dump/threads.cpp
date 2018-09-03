#include <pthread.h>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[]);	//arguments to main will be size of Lattice
/*
In main function we define a temp set and create as many threads as # of temp points, each seeded with a different random number.
Every thread generates a random lattice that updates all spins in one "step"
*/