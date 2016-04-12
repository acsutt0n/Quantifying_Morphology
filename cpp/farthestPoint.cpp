// Compute the farthest possible distance between the given points
//


#include <eigen3/Eigen/Dense>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>


using namespace Eigen;
using namespace std;


template<std::size_t n>
void printBoolSize(bool (&boolarray)[n]){
  cout << "size of array is: " << sizeof(boolarray) << endl;
  cout << "size of a given element is: " << sizeof(boolarray[0]) << endl;
}


int main() {
  int myinput = 1000;
  /*
  string input = "";
  cout << "how many elements?  " << endl;
  getline(cin, input);
  stringstream myStream(input);
  myStream >> myinput; 
  bool betch[myinput] = {false};
  for (int i = 0; i < myinput; i++) {
    betch.push_back(33);
  } */
  bool betch[10000][50] = {false};
  printBoolSize(betch);
  return 0;
}

