#include "spherical_harmonics.hpp"

#include <complex>
#include <iostream>

typedef double sh_digit;

namespace math {
namespace {
  // Hardcoded spherical harmonic functions for low orders (l is first number
  // and m is second number (sign encoded as preceeding 'p' or 'n')).
  //
  // As polynomials they are evaluated more efficiently in cartesian coordinates,
  // assuming that @d is unit. This is not verified for efficiency.

  // https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
  std::complex<sh_digit> HardcodedSH6m6(const dir &d) {
    // N[1/64 *Sqrt[3003/Pi], 16] e-ipsi^6 sinTheta^6
    sh_digit tmp = std::pow(d.sinTheta, 6);
    std::complex<sh_digit> emipsi(d.cosPsi, -d.sinPsi);
    std::complex<sh_digit> ans(std::pow(emipsi, 6));
    ans *= tmp * 0.4830841135800662;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH6m5(const dir &d) {
    // N[3/32 Sqrt[1001/Pi], 16] e-ipsi^5 sinTheta^5 cosTheta
    sh_digit tmp = std::pow(d.sinTheta, 5) * d.cosTheta;
    std::complex<sh_digit> emipsi(d.cosPsi, -d.sinPsi);
    std::complex<sh_digit> ans(std::pow(emipsi, 5));
    ans *= tmp * 1.673452458100098;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH6m4(const dir &d) {
    // N[3/32 Sqrt[91/Pi/2], 16] e-ipsi^4 sinTheta^4 (11 cosTheta^2-1)
    sh_digit tmp = std::pow(d.sinTheta, 4) * (11 * d.cosTheta * d.cosTheta - 1);
    std::complex<sh_digit> emipsi(d.cosPsi, -d.sinPsi);
    std::complex<sh_digit> ans(std::pow(emipsi, 4));
    ans *= tmp * 0.3567812628539980;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH6m3(const dir &d) {
    // N[1/32 Sqrt[1365/Pi], 16] e-ipsi^3 sinTheta^3 (11 cosTheta^3-3 cosTheta)
    sh_digit tmp = std::pow(d.sinTheta, 3) * (11 * std::pow(d.cosTheta, 3) - 3 * d.cosTheta);
    std::complex<sh_digit> emipsi(d.cosPsi, -d.sinPsi);
    std::complex<sh_digit> ans(std::pow(emipsi, 3));
    ans *= tmp * 0.6513904858677157;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH6m2(const dir &d) {
    // N[1/64 Sqrt[1365/Pi], 16] e-ipsi^2 sinTheta^2 (33 cosTheta^4 -18 cosTheta^2 +1)
    sh_digit tmp = std::pow(d.sinTheta, 2) * (33 * std::pow(d.cosTheta, 4) - 18 * d.cosTheta * d.cosTheta + 1);
    std::complex<sh_digit> emipsi(d.cosPsi, -d.sinPsi);
    std::complex<sh_digit> ans(std::pow(emipsi, 2));
    ans *= tmp * 0.3256952429338579;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH6m1(const dir &d) {
    // N[1/16 Sqrt[273/Pi/2], 16] e-ipsi^1 sinTheta^1 (33 cosTheta^5 -30 cosTheta^3 +5 cosTheta)
    sh_digit tmp = d.sinTheta * (33 * std::pow(d.cosTheta, 5) - 30 * std::pow(d.cosTheta, 3) + 5 * d.cosTheta);
    std::complex<sh_digit> ans(d.cosPsi, -d.sinPsi);
    ans *= tmp * 0.4119755163011408;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH60(const dir &d) {
    // N[1/32 Sqrt[13/Pi], 16] (231 cosTheta^6 -315 cosTheta^4 +105 cosTheta^2 -5)
    sh_digit tmp = (231 * std::pow(d.cosTheta, 6) - 315 * std::pow(d.cosTheta, 4) + 105 * std::pow(d.cosTheta, 2) - 5);
    std::complex<sh_digit> ans(1, 0);
    ans *= tmp * 0.06356920226762843;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH61(const dir &d) {
    // N[-1/16 Sqrt[273/Pi/2], 16] eipsi^1 sinTheta (33 cosTheta^5 -30 cosTheta^3 +5 cosTheta)
    sh_digit tmp = d.sinTheta * (33 * std::pow(d.cosTheta, 5) - 30 * std::pow(d.cosTheta, 3) + 5 * d.cosTheta);
    std::complex<sh_digit> ans(d.cosPsi, d.sinPsi);
    ans *= tmp * -0.4119755163011408;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH62(const dir &d) {
    // N[1/64 Sqrt[1365/Pi], 16] eipsi^2 sinTheta^2 (33 cosTheta^4 -18 cosTheta^2 +1)
    sh_digit tmp = std::pow(d.sinTheta, 2) * (33 * std::pow(d.cosTheta, 4) - 18 * d.cosTheta * d.cosTheta + 1);
    std::complex<sh_digit> eipsi(d.cosPsi, d.sinPsi);
    std::complex<sh_digit> ans(std::pow(eipsi, 2));
    ans *= tmp * 0.3256952429338579;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH63(const dir &d) {
    // N[-1/32 Sqrt[1365/Pi], 16] eipsi^3 sinTheta^3 (11 cosTheta^3-3 cosTheta)
    sh_digit tmp = std::pow(d.sinTheta, 3) * (11 * std::pow(d.cosTheta, 3) - 3 * d.cosTheta);
    std::complex<sh_digit> eipsi(d.cosPsi, d.sinPsi);
    std::complex<sh_digit> ans(std::pow(eipsi, 3));
    ans *= tmp * -0.6513904858677157;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH64(const dir &d) {
    // N[3/32 Sqrt[91/Pi/2], 16] eipsi^4 sinTheta^4 (11 cosTheta^2-1)
    sh_digit tmp = std::pow(d.sinTheta, 4) * (11 * d.cosTheta * d.cosTheta - 1);
    std::complex<sh_digit> eipsi(d.cosPsi, d.sinPsi);
    std::complex<sh_digit> ans(std::pow(eipsi, 4));
    ans *= tmp * 0.3567812628539980;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH65(const dir &d) {
    // N[-3/32 Sqrt[1001/Pi], 16] eipsi^5 sinTheta^5 cosTheta
    sh_digit tmp = std::pow(d.sinTheta, 5) * d.cosTheta;
    std::complex<sh_digit> eipsi(d.cosPsi, d.sinPsi);
    std::complex<sh_digit> ans(std::pow(eipsi, 5));
    ans *= tmp * -1.673452458100098;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH66(const dir &d) {
    // N[1/64 *Sqrt[3003/Pi], 16] eipsi^6 sinTheta^6
    sh_digit tmp = std::pow(d.sinTheta, 6);
    std::complex<sh_digit> eipsi(d.cosPsi, d.sinPsi);
    std::complex<sh_digit> ans(std::pow(eipsi, 6));
    ans *= tmp * 0.4830841135800662;
    return ans;
  }

  // start Y4m part
  std::complex<sh_digit> HardcodedSH4m4(const dir &d) {
    // N[3/16 Sqrt[35/Pi/2], 16] e-ipsi^4 sinTheta^4
    sh_digit tmp = std::pow(d.sinTheta, 4);
    std::complex<sh_digit> emipsi(d.cosPsi, -d.sinPsi);
    std::complex<sh_digit> ans(std::pow(emipsi, 4));
    ans *= tmp * 0.4425326924449826;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH4m3(const dir &d) {
    // N[3/8 Sqrt[35/Pi], 16] e-ipsi^3 sinTheta^3 cosTheta
    sh_digit tmp = std::pow(d.sinTheta, 3) * d.cosTheta;
    std::complex<sh_digit> emipsi(d.cosPsi, -d.sinPsi);
    std::complex<sh_digit> ans(std::pow(emipsi, 3));
    ans *= tmp * 1.251671470898352;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH4m2(const dir &d) {
    // N[3/8 Sqrt[5/Pi/2], 16] e-ipsi^2 sinTheta^2 (7 cosTheta^2 -1)
    sh_digit tmp = std::pow(d.sinTheta, 2) * (7 * std::pow(d.cosTheta, 2) - 1);
    std::complex<sh_digit> emipsi(d.cosPsi, -d.sinPsi);
    std::complex<sh_digit> ans(std::pow(emipsi, 2));
    ans *= tmp * 0.3345232717786446;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH4m1(const dir &d) {
    // N[3/8 Sqrt[5/Pi], 16] e-ipsi^1 sinTheta^1 (7 cosTheta^3 -3 cosTheta)
    sh_digit tmp = d.sinTheta * (7 * std::pow(d.cosTheta, 3) - 3 * d.cosTheta);
    std::complex<sh_digit> ans(d.cosPsi, -d.sinPsi);
    ans *= tmp * 0.4730873478787800;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH40(const dir &d) {
    // N[3/16 Sqrt[1/Pi], 16] (35 cosTheta^4 -30 cosTheta^2 +3)
    sh_digit tmp = (35 * std::pow(d.cosTheta, 4) - 30 * std::pow(d.cosTheta, 2) + 3);
    std::complex<sh_digit> ans(1, 0);
    ans *= tmp * 0.1057855469152043;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH41(const dir &d) {
    // N[-3/8 Sqrt[5/Pi], 16] eipsi^1 sinTheta^1 (7 cosTheta^3 -3 cosTheta)
    sh_digit tmp = d.sinTheta * (7 * std::pow(d.cosTheta, 3) - 3 * d.cosTheta);
    std::complex<sh_digit> ans(d.cosPsi, d.sinPsi);
    ans *= tmp * -0.4730873478787800;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH42(const dir &d) {
    // N[3/8 Sqrt[5/Pi/2], 16] eipsi^2 sinTheta^2 (7 cosTheta^2 -1)
    sh_digit tmp = std::pow(d.sinTheta, 2) * (7 * std::pow(d.cosTheta, 2) - 1);
    std::complex<sh_digit> eipsi(d.cosPsi, d.sinPsi);
    std::complex<sh_digit> ans(std::pow(eipsi, 2));
    ans *= tmp * 0.3345232717786446;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH43(const dir &d) {
    // N[-3/8 Sqrt[35/Pi], 16] eipsi^3 sinTheta^3 cosTheta
    sh_digit tmp = std::pow(d.sinTheta, 3) * d.cosTheta;
    std::complex<sh_digit> eipsi(d.cosPsi, d.sinPsi);
    std::complex<sh_digit> ans(std::pow(eipsi, 3));
    ans *= tmp * -1.251671470898352;
    return ans;
  }
  std::complex<sh_digit> HardcodedSH44(const dir &d) {
    // N[3/16 Sqrt[35/Pi/2], 16] e-ipsi^4 sinTheta^4
    sh_digit tmp = std::pow(d.sinTheta, 4);
    std::complex<sh_digit> eipsi(d.cosPsi, d.sinPsi);
    std::complex<sh_digit> ans(std::pow(eipsi, 4));
    ans *= tmp * 0.4425326924449826;
    return ans;
  }

}  // namespace

std::complex<sh_digit> Y6m(int m, const dir &d) {
  switch (m) {
    case -6:
      return HardcodedSH6m6(d);
    case -5:
      return HardcodedSH6m5(d);
    case -4:
      return HardcodedSH6m4(d);
    case -3:
      return HardcodedSH6m3(d);
    case -2:
      return HardcodedSH6m2(d);
    case -1:
      return HardcodedSH6m1(d);
    case 0:
      return HardcodedSH60(d);
    case 1:
      return HardcodedSH61(d);
    case 2:
      return HardcodedSH62(d);
    case 3:
      return HardcodedSH63(d);
    case 4:
      return HardcodedSH64(d);
    case 5:
      return HardcodedSH65(d);
    case 6:
      return HardcodedSH66(d);
  }

  // This is unreachable given the CHECK's above but the compiler can't tell.
  return 0.0;
}

sh_digit Y6m_sqnorm(int m, const dir &d) {
  auto tmp = Y6m(m, d);
  return tmp.real() * tmp.real() + tmp.imag() * tmp.imag();
}

std::complex<sh_digit> Y4m(int m, const dir &d) {
  switch (m) {
    case -4:
      return HardcodedSH4m4(d);
    case -3:
      return HardcodedSH4m3(d);
    case -2:
      return HardcodedSH4m2(d);
    case -1:
      return HardcodedSH4m1(d);
    case 0:
      return HardcodedSH40(d);
    case 1:
      return HardcodedSH41(d);
    case 2:
      return HardcodedSH42(d);
    case 3:
      return HardcodedSH43(d);
    case 4:
      return HardcodedSH44(d);
  }

  // This is unreachable given the CHECK's above but the compiler can't tell.
  return 0.0;
}

sh_digit Y4m_sqnorm(int m, const dir &d) {
  auto tmp = Y4m(m, d);
  return tmp.real() * tmp.real() + tmp.imag() * tmp.imag();
}
}  // namespace math

// In[8]:= N[1/64 *Sqrt[3003/Pi], 16]

// Out[8]= 0.4830841135800662

// In[9]:= N[3/32 Sqrt[1001/Pi], 16]

// Out[9]= 1.673452458100098

// In[17]:= N[3/32 Sqrt[91/Pi/2], 16]

// Out[17]= 0.3567812628539980

// In[11]:= N[1/32 Sqrt[1365/Pi], 16]

// Out[11]= 0.6513904858677157

// In[12]:= N[1/64 Sqrt[1365/Pi], 16]

// Out[12]= 0.3256952429338579

// In[13]:= N[1/16 Sqrt[273/Pi/2], 16]

// Out[13]= 0.4119755163011408

// In[14]:= N[1/32 Sqrt[13/Pi], 16]

// Out[14]= 0.06356920226762843

// In[6]:= N[3/16 Sqrt[35/Pi/2], 16]

// Out[6]= 0.4425326924449826

// In[7]:= N[3/8 Sqrt[35/Pi], 16]

// Out[7]= 1.251671470898352

// In[8]:= N[3/8 Sqrt[5/Pi/2], 16]

// Out[8]= 0.3345232717786446

// In[9]:= N[3/8 Sqrt[5/Pi], 16]

// Out[9]= 0.4730873478787800

// In[10]:= N[3/16 Sqrt[1/Pi], 16]

// Out[10]= 0.1057855469152043