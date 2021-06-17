#ifndef WIGNER3J_WIGNER3J_HPP
#define WIGNER3J_WIGNER3J_HPP

// Calculate the wigner3j symbol for Wl.
// When 3js are equal, the coffs hold an analytic form.
// The analytic form is obtained by using Mathematica
// Assuming[Element[a | b, Integers], ThreeJSymbol[{6, a}, {6, b}, {6, -a - b}]]
// and
// Assuming[Element[a | b, Integers], ThreeJSymbol[{4, a}, {4, b}, {4, -a - b}]]

#include <cmath>
#include <vector>

namespace math {

const int kCacheSize = 16;
double Factorial(int x);

// For l=6, the range of m_i is l1 in [-6,6] and l2 in [-6,6] and l1+l2 in [-6,6]
double w3_6(int a, int b);

double w3_4(int a, int b);
}  // namespace math

#endif  //WIGNER3J_WIGNER3J_HPP
