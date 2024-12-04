
// simple demo program for svd3x3

#include <iostream>
#include "svd3x3.h"

#include <chrono>

int main() {
  //double mat[9] = {-1, 4, -7, 4, 4, 4, -2, 2, 7};
  double mat[9] = {8, 4, 17, 3, 9, 2, 27, 2, 1};

  std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
  start = std::chrono::high_resolution_clock::now();

  double U[9], V[9], s[3];
  svd3x3::decomp<double> svd3;
  svd3(mat, U, s, V);

  stop = std::chrono::high_resolution_clock::now();
  auto micros = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Time taken by svd3x3: " << micros.count() << " microseconds" << std::endl << std::endl;

  printf("U\n");
  printf("%g %g %g\n", U[0], U[1], U[2]);
  printf("%g %g %g\n", U[3], U[4], U[5]);
  printf("%g %g %g\n", U[6], U[7], U[8]);
  printf("s: %g %g %g\n", s[0], s[1], s[2]);
  printf("V\n");
  printf("%g %g %g\n", V[0], V[1], V[2]);
  printf("%g %g %g\n", V[3], V[4], V[5]);
  printf("%g %g %g\n", V[6], V[7], V[8]);

  // check
  svd3.compose(U, s, V, mat);
  printf("Check U*S*Vt:\n");
  printf("%g %g %g\n", mat[0], mat[1], mat[2]);
  printf("%g %g %g\n", mat[3], mat[4], mat[5]);
  printf("%g %g %g\n", mat[6], mat[7], mat[8]);
}
