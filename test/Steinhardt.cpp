#include "Steinhardt.hpp"

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <voro++/voro++.hh>

#include "Frame.hpp"
#include "spherical_harmonics.hpp"
#include "wigner3j.hpp"

#define Debugql 0

const double PI = std::acos(-1);

// should move to  util
inline double periodic_distance(double x, double lx)
{
  return (std::abs(x) > 0.5 * lx) ? ((x > 0) ? x - lx : x + lx) : x;
}

// return the Y_l-m part
// Only positive m is calculated, here return the minus part
std::complex<double> q_lm(const std::vector<std::complex<double>> &qllist, const int &m)
{
  if (m < 0)
    // -m to keep it inside
    return std::pow(-1, -m % 2) * std::conj(qllist[-m]);
  else
    return qllist[m];
}
Table qw(const Frame &data)
{
  // calculated properties
  std::vector<dtype> q6(data.particleN), q4(data.particleN), w6(data.particleN), w4(data.particleN);

  double Cn = std::cbrt(data.particleN / 8 / data.xl / data.yl / data.zl);
  voro::container con(data.boxXL, data.boxXH, data.boxYL, data.boxYH, data.boxZL, data.boxZH,
                      int(std::ceil(Cn * data.xl)), int(std::ceil(Cn * data.yl)), int(std::ceil(Cn * data.zl)),
                      true, true, true, 8);
  // id start from 1, so here should minus 1 to keep size tight
  for (size_t i = 0; i < data.particleN; i++)
    con.put(i, data.x[i], data.y[i], data.z[i]);
  // con.put(id[i] - 1, x[i], y[i], z[i]);

  std::vector<std::vector<std::complex<dtype>>> q6List(data.particleN, std::vector<std::complex<dtype>>(7));
  std::vector<std::vector<std::complex<dtype>>> q4List(data.particleN, std::vector<std::complex<dtype>>(5));

#if Debugql
  if (q6List[0].size() != 7)
    throw std::runtime_error("7");
#endif
  voro::voronoicell_neighbor cell_neighbor;
  std::vector<int> neighbor;
  // store total area for each particle
  std::vector<double> totalArea(data.particleN, 0);
  std::vector<double> cellArea;

#if Debugql
  std::vector<double> totalAreaDebug(data.particleN, 0);

  std::cout << "begin cloop\n";
  bool once = true;
  size_t N = 0;
#endif

  voro::c_loop_all cloop(con);
  if (cloop.start())
    do
      // compute_cell input voronoicell_base
      //              voronoicell_neighbor inheritate from that
      // return info to cell_neighbor
      if (con.compute_cell(cell_neighbor, cloop))
      {
        // the current  id of particle
        int id = cloop.pid();

        // Gather information about the computed Voronoi cell
        // return neighbor info of one particle to a vector
        cell_neighbor.neighbors(neighbor);
        cell_neighbor.face_areas(cellArea);

#if Debugql
        for (size_t i = 0; i < cellArea.size(); i++)
          totalAreaDebug[id] += cellArea[i];
        if (once || N == 179)
          std::cout << "cell area right\n";
#endif

        // Loop over all faces of the Voronoi cell
        for (size_t i = 0; i < neighbor.size(); i++)
        {
          // Skip if the neighbor information is smaller than
          // this particle's ID, to avoid double counting.
          if (neighbor[i] > id)
          {
            // add area to total area
            totalArea[id] += cellArea[i];
            totalArea[neighbor[i]] += cellArea[i];
            // consider the periodic boundary condition
            math::dir v(periodic_distance(data.x[id] - data.x[neighbor[i]], data.xl),
                        periodic_distance(data.y[id] - data.y[neighbor[i]], data.yl),
                        periodic_distance(data.z[id] - data.z[neighbor[i]], data.zl));

#if Debugql
            if (once || N == 179)
              std::cout << "dir init right\n";
#endif

            std::complex<dtype> tmp = cellArea[i] * math::Y6m(0, v);
            q6List[id][0] += tmp;
            // q6List.at(id).at(0) += tmp;
            // q6List.at(neighbor.at(i)).at(0) += tmp;
            q6List[neighbor[i]][0] += tmp;
            tmp = cellArea[i] * math::Y4m(0, v);
            q4List[id][0] += tmp;
            q4List[neighbor[i]][0] += tmp;

#if Debugql
            if (once || N == 179)
              std::cout << "Ylm added right\n";
#endif
            for (int m = 1; m < 7; m++)
            {
              std::complex<dtype> tmp = cellArea[i] * math::Y6m(m, v);
              q6List[id][m] += tmp;
              q6List[neighbor[i]][m] += tmp;
#if Debugql
              if (once)
                std::cout << "Y6m right\n";
#endif
            }
            for (int m = 1; m < 5; m++)
            {
              std::complex<dtype> tmp = cellArea[i] * math::Y4m(m, v);
              q4List[id][m] += tmp;
              q4List[neighbor[i]][m] += tmp;
            }
          }
#if Debugql
          if (once || N == 179)
            std::cout << "one loop complete\n";
#endif
          // Skip to the next entry in the face vertex list
          // j += f_vert[j] + 1;
        }
#if Debugql
        if (once || N == 179)
          std::cout << "one cell complete\n";
#endif
#if Debugql
        once = false;
        N += 1;
        std::cout << N << "\n";
#endif
      }
    while (cloop.inc());
#if Debugql
  std::cout << "compute complete\n";
#endif
#if Debugql
  for (size_t i = 0; i < totalArea.size(); i++)
    if (std::abs(totalAreaDebug[i] - totalArea[i]) > 1e-8)
    {
      std::cout << "area is not match" << std::endl;
      std::cout << totalAreaDebug[i] << " " << totalArea[i] << std::endl;
      break;
    };

  std::cout << "area check complete\n";
#endif

  const dtype constCoffQ6 = std::sqrt(4 * PI / (2 * 6 + 1));
  const dtype constCoffQ4 = std::sqrt(4 * PI / (2 * 4 + 1));

  for (size_t i = 0; i < data.particleN; i++)
  {
    q6[i] += q6List[i][0].real() * q6List[i][0].real() + q6List[i][0].imag() * q6List[i][0].imag();
    for (int m = 1; m < 7; m++)
      q6[i] += 2 * (q6List[i][m].real() * q6List[i][m].real() + q6List[i][m].imag() * q6List[i][m].imag());
    q6[i] = std::sqrt(q6[i]);
    q6[i] *= constCoffQ6 / totalArea[i];
  }
  for (size_t i = 0; i < data.particleN; i++)
  {
    q4[i] += q4List[i][0].real() * q4List[i][0].real() + q4List[i][0].imag() * q4List[i][0].imag();
    for (int m = 1; m < 5; m++)
      q4[i] += 2 * (q4List[i][m].real() * q4List[i][m].real() + q4List[i][m].imag() * q4List[i][m].imag());
    q4[i] = std::sqrt(q4[i]);
    q4[i] *= constCoffQ4 / totalArea[i];
  }
  for (size_t i = 0; i < data.particleN; i++)
  {
    for (int a = -6; a < 7; a++)
    {
      for (int b = -6; b < 7; b++)
      {
        if (std::abs(a + b) < 7)
          w6[i] += math::w3_6(a, b) * (q_lm(q6List[i], a) * q_lm(q6List[i], b) * q_lm(q6List[i], -a - b)).real();
      }
    }
    w6[i] /= std::pow(q6[i], 3) * std::pow(totalArea[i], 3);
    w6[i] *= std::pow(constCoffQ6, 3);
    for (int a = -4; a < 5; a++)
    {
      for (int b = -4; b < 5; b++)
      {
        if (std::abs(a + b) < 5)
          w4[i] += math::w3_4(a, b) * (q_lm(q4List[i], a) * q_lm(q4List[i], b) * q_lm(q4List[i], -a - b)).real();
      }
    }
    w4[i] /= std::pow(q4[i], 3) * std::pow(totalArea[i], 3);
    w4[i] *= std::pow(constCoffQ4, 3);
  }

  return Table({q6, q4, w6, w4});
}
