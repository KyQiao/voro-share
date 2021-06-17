// Calculate the wigner3j symbol for Wl.
// When 3js are equal, the coffs hold an analytic form.
// The analytic form is obtained by using Mathematica
// Assuming[Element[a | b, Integers], ThreeJSymbol[{6, a}, {6, b}, {6, -a - b}]]
// and
// Assuming[Element[a | b, Integers], ThreeJSymbol[{4, a}, {4, b}, {4, -a - b}]]

#include "wigner3j.hpp"

#include <cmath>
#include <vector>

namespace math {
double Factorial(int x) {
  const double factorial_cache[kCacheSize] = {1, 1, 2, 6, 24, 120, 720, 5040,
                                              40320, 362880, 3628800, 39916800,
                                              479001600, 6227020800,
                                              87178291200, 1307674368000};
  return factorial_cache[x];
}

// For l=6, the range of m_i is l1 in [-6,6] and l2 in [-6,6] and l1+l2 in [-6,6]
double w3_6(int a, int b) {
  // (1/(2 Sqrt[46189]))((-1)^-a)
  // (-28800+15568 a^2
  // -1190 a^4+22 a^6
  // +58520 a b
  // -8113 a^3 b
  // +209 a^5 b
  // +58520 b^2
  // -21679 a^2 b^2
  // +855 a^4 b^2
  // -27132 a b^3
  // +1938 a^3 b^3
  // -13566 b^4
  // +2584 a^2 b^4
  // +1938 a b^5
  // +646 b^6)
  // +sqrt(((6-a)! (6+a)!)/((6-b)! (6-a-b)! (6+b)! (6+a+b)!))

  double ans = 0.002326487505282269;
  int c = -a - b;
  if ((std::abs(a) > 6) || (std::abs(b) > 6) || (std::abs(c) > 6))
    return 0;
  else {
    double tmp = -28800 + 15568 * std::pow(a, 2)
                 - 1190 * std::pow(a, 4)
                 + 22 * std::pow(a, 6)
                 + 58520 * a * b
                 - 8113 * std::pow(a, 3) * b
                 + 209 * std::pow(a, 5) * b
                 + 58520 * std::pow(b, 2)
                 - 21679 * std::pow(a, 2) * std::pow(b, 2)
                 + 855 * std::pow(a, 4) * std::pow(b, 2)
                 - 27132 * a * std::pow(b, 3)
                 + 1938 * std::pow(a, 3) * std::pow(b, 3)
                 - 13566 * std::pow(b, 4)
                 + 2584 * std::pow(a, 2) * std::pow(b, 4)
                 + 1938 * a * std::pow(b, 5)
                 + 646 * std::pow(b, 6);
    double tmp2 = Factorial(6 - a) * Factorial(6 + a) / (Factorial(6 - b) * Factorial(6 + c) * Factorial(6 + b) * Factorial(6 - c));
    tmp2 = std::sqrt(tmp2);
    ans *= tmp * tmp2;
    return (a % 2 == 0) ? ans : -ans;
  }
}

double w3_4(int a, int b) {
  // (1/(3 Sqrt[2002]))((-1)^-a)
  // (432-230 a^2
  // +14 a^4
  // -845 a b
  // +91 a^3 b
  // -845 b^2
  // +234 a^2 b^2
  // +286 a b^3
  // +143 b^4)
  double ans = 0.007449835937794570;

  int c = -a - b;
  if ((std::abs(a) > 4) || (std::abs(b) > 4) || (std::abs(c) > 4))
    return 0;
  else {
    double tmp = 432 - 230 * std::pow(a, 2)
                 + 14 * std::pow(a, 4)
                 - 845 * a * b
                 + 91 * std::pow(a, 3) * b
                 - 845 * std::pow(b, 2)
                 + 234 * std::pow(a, 2) * std::pow(b, 2)
                 + 286 * a * std::pow(b, 3)
                 + 143 * std::pow(b, 4);
    double tmp2 = Factorial(4 - a) * Factorial(4 + a) / (Factorial(4 - b) * Factorial(4 + c) * Factorial(4 + b) * Factorial(4 - c));
    tmp2 = std::sqrt(tmp2);
    ans *= tmp * tmp2;
    return (a % 2 == 0) ? ans : -ans;
  }
}

}  // namespace math
