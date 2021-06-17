#ifndef SH_SPHERICAL_HARMONICS_H_
#define SH_SPHERICAL_HARMONICS_H_

#include <cmath>
#include <complex>

typedef double sh_digit;

namespace math {

// A unit vector, represent the angle theta and psi.
struct dir {
  // double x,y,z;
  sh_digit sinTheta, cosTheta, sinPsi, cosPsi;
  std::complex<sh_digit> eipsi;
  dir(sh_digit x, sh_digit y, sh_digit z) {
    sh_digit rxy2 = x * x + y * y;
    sh_digit rxy = std::sqrt(rxy2);
    sh_digit r2 = rxy2 + z * z;
    sh_digit r = std::sqrt(r2);
    sinTheta = rxy / r;
    cosTheta = z / r;
    if (rxy == 0) {
      // deal with the nan situation
      // when r= 0 0 1, sh should return 0;
      sinPsi = 0;
      cosPsi = 0;
    } else {
      sinPsi = x / rxy;
      cosPsi = y / rxy;
    }
  }
};

// Evaluate the spherical harmonic basis function of degree @l and order @m
// for the given direction vector, @dir.
// For low values of @l this will use a hard-coded function, otherwise it
// will fallback to EvalSHSlow that uses a recurrence relation to support all l.
std::complex<sh_digit> Y6m(int m, const dir &d);

std::complex<sh_digit> Y4m(int m, const dir &d);

// return the squared norm of Y6m
sh_digit Y6m_sqnorm(int m, const dir &d);

// return the squared norm of Y4m
sh_digit Y4m_sqnorm(int m, const dir &d);
}  // namespace sh

#endif  // SH_SPHERICAL_HARMONICS_H_
