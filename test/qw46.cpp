#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include "Frame.hpp"
#include "Steinhardt.hpp"

double mean(std::vector<double> v)
{
  size_t size = v.size();
  double ans = 0;
  for (size_t i = 0; i < size; i++)
    ans += v[i];
  ans /= size;
  return ans;
}

int main(int argc, char const *argv[])
{
  if ((argc != 2))
  {
    std::cout << "input number should be one" << std::endl;
    std::cout << "qw fcc.xyz\n"
              << std::endl;
    return 0;
  }

  std::string filename = argv[1];

  auto data = Frame();
  data.read(filename);
  // q6 q4 w6 w4
  Table result = qw(data);

  if (filename == "fcc")
  {
    double q6 = mean(result.data[0]), q4 = mean(result.data[1]), w6 = mean(result.data[2]), w4 = mean(result.data[3]);
    if ((q6 != 0.574524) && (q4 != 0.190941) && (w6 != -0.0131606) && (w4 != -0.159317))
      return -1;
  }

  return 0;
}
