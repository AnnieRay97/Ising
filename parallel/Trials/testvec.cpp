#include <iostream>
#include <vector>
#include <array>
using namespace std;

int main()
{
	//vector <int> lattices;
	// for (int i=0; i< 5; i++)
	// 	lattices.push_back(i);
	// for (int i=0; i<5; i++)
	// 	cout << lattices.at(i);
	cout << "input size" << endl;
	int b;
	cin >> b;
	array<int, b> a = {1,2,3};
	cout << a[2];
	return 0;
}