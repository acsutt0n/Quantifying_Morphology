// Represent a full-sized matrix

#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <eigen3/Eigen/Dense>
#include <ctime>
#include <cstdlib>


using namespace std;
using namespace Eigen;



MatrixXd getMatrix(int cols=1000, int rows=2000) {
  // int const RAND_MAX = 2;
  srand(time(0)); // seed rand num gen with current time
  MatrixXd mat(rows, cols);
  cout << "size of newly-allocated matrix: " << mat.size() << endl;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      mat(i,j) = rand(); // random value on [0:RAND_MAX]
    }
  }
  cout << "size of matrix after populating: " << mat.size() << endl;
  return mat;
}


int main() {
  // boost::dynamic_bitset<>
  
  MatrixXd m = getMatrix(2000,2000);
  return 0;
}

