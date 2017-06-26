#include <iostream>

using namespace std;

extern "C"
void func()
{
  cout << "func called!" << endl;
}

int a = 5;

