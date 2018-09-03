#include <iostream>
#include <complex>
#include <cmath>
#define L 5
using namespace std;
typedef complex<double> dcomp;

main() {
  dcomp i, a;
  int N = 50;
  double pi;
  pi = 2 * asin(1);
  i = -1;
  i = sqrt(i);
  a = exp(2*(pi/N)*i);
  cout << a << endl;
  cout << abs(a) << endl;
  a += exp (2*(pi/N)*i);
  cout << a << endl;
  a /= L;
  cout << a << endl;
  cout << abs (a);
  return 0;
} 
