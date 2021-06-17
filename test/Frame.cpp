#include "Frame.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "type.hpp"
#include "zstr.hpp"

#define DebugInput 0

// https://codereview.stackexchange.com/questions/206686/removing-by-indices-several-elements-from-a-vector
template <typename T, typename Iter>
void removeIndicesFromVector(std::vector<T>& v, Iter begin, Iter end)
// requires std::is_convertible_v<std::iterator_traits<Iter>::value_type, std::size_t>
{
  assert(std::is_sorted(begin, end));
  auto rm_iter = begin;
  std::size_t current_index = 0;

  const auto pred = [&](const T&) {
    // any more to remove?
    if (rm_iter == end) {
      return false;
    }
    // is this one specified?
    if (*rm_iter == current_index++) {
      return ++rm_iter, true;
    }
    return false;
  };

  v.erase(std::remove_if(v.begin(), v.end(), pred), v.end());
}

template <typename T, typename S>
// requires std::is_convertible_v<S::value_type, std::size_t>
void removeIndicesFromVector(std::vector<T>& v, const S& rm) {
  using std::begin;
  using std::end;
  assert(std::is_sorted(begin(rm), end(rm)));
  return removeIndicesFromVector(v, begin(rm), end(rm));
}

template <typename fileInput>
void Frame::readLammpsBase(fileInput& file, bool readAttr) {
  // a temp string for storing lines
  std::string line;
  std::stringstream ss;

  // check the first line
  getline(file, line);
  if (!(line.substr(0, 4) == "ITEM")) {
    std::cout << "start with " << line << std::endl;
    throw std::runtime_error("this is not a LAMMPS dump file!");
  }

  // skip lines in file
  auto skipLines = [&file, &line](int i) {
    while (i-- > 0)
      getline(file, line);
  };

  // get string into ss
  auto ssGet = [&file, &line, &ss]() {
    getline(file, line);
    // https://stackoverflow.com/questions/20731/how-do-you-clear-a-stringstream-variable
    // https://stackoverflow.com/questions/13891856/cannot-overwrite-stringstream-variable-with-a-new-value
    ss.str("");
    ss.clear();
    ss << line;
    // using sstream is quite tricky. above is the way to rewrite buffer
  };

  ssGet();
  ss >> this->timestep;
  // this->timestep = std::stoul(ss.str());

  // skip particle line
  skipLines(1);
  ssGet();
  ss >> this->particleN;
  // this->particleN = std::stoul(ss.str());

  // skip box line
  skipLines(1);

  ssGet();
  ss >> this->boxXL >> this->boxXH;
  ssGet();
  ss >> this->boxYL >> this->boxYH;
  ssGet();
  ss >> this->boxZL >> this->boxZH;

  // parse input order
  getline(file, line);
  // the first 11 char is ITEM: ATOMS
  ss.str("");
  ss.clear();
  ss << line.substr(12);
#if DebugInput
  std::cout << ss.str() << std::endl;
#endif
  std::vector<std::function<void(size_t, std::stringstream&)>> callOrder;
  callOrder.reserve(20);
  // index for dict lookup
  size_t index = 0;

  // https://stackoverflow.com/questions/8046357/how-do-i-check-if-a-stringstream-variable-is-empty-null
  std::string tmp;
  // while (ss.rdbuf()->in_avail()!=0) {
  while (std::getline(ss, tmp, ' ')) {
#if DebugInput
    std::cout << "tmp is " << tmp << std::endl;
    std::cout << "and ss is " << ss.str() << std::endl;
#endif
    if (tmp == "x") {
      this->x.resize(this->particleN);
      callOrder.push_back([this](size_t index, std::stringstream& s) { s >> this->x[index]; });
    } else if (tmp == "y") {
      this->y.resize(this->particleN);
      callOrder.push_back([this](size_t index, std::stringstream& s) { s >> this->y[index]; });
    } else if (tmp == "z") {
      this->z.resize(this->particleN);
      callOrder.push_back([this](size_t index, std::stringstream& s) { s >> this->z[index]; });
    } else if (tmp == "id") {
      this->id.resize(this->particleN);
      callOrder.push_back([this](size_t index, std::stringstream& s) { s >> this->id[index]; });
    } else if (tmp == "type") {
      this->type.resize(this->particleN);
      callOrder.push_back([this](size_t index, std::stringstream& s) { s >> this->type[index]; });
    } else {
      if (readAttr) {
        size_t col = (this->attr_table).cols();
#if DebugInput
        std::cout << "col is " << col << std::endl;
#endif
        // (this->attr_table).data.push_back({});
        // (this->attr_table).setcols(5);
        (this->attr_table).setcols(col + 1);
#if DebugInput
        std::cout << "binding cols " << std::endl;
#endif
        callOrder.push_back([col, this](size_t index, std::stringstream& s) { s >> (this->attr_table).data[col][index]; });
        this->attr_order += tmp + " ";
        this->attr_index[tmp] = index++;
      }
    }
  }

#if DebugInput
  std::cout << "function vector loaded " << std::endl;
#endif
  // resize space for data for faster load

  (this->attr_table).setrows(this->particleN);
#if DebugInput
  std::cout << "attribute table seted " << std::endl;
#endif
  index = 0;
  while (getline(file, line)) {
    std::stringstream tmp(line);
#if DebugInput
    std::cout << tmp.str() << std::endl;
#endif
    for (auto& f : callOrder) {
      f(index, tmp);
    }
    index++;
  }
}

// read
void Frame::readLammps(const std::string& fileName, compressType compress, bool readAttr) {
  // store read order from file
  std::vector<double> readOrder{};

  if (compress == gz) {
    // read binary file into gz
    // must use ifstream here, other port won't be matched
    zstr::ifstream file(fileName, std::ios_base::in | std::ios_base::binary);
    readLammpsBase<zstr::ifstream>(file, readAttr);

  } else if (compress == none) {
    // read plain file
    std::ifstream file(fileName, std::ios_base::in);
    readLammpsBase<decltype(file)>(file, readAttr);
    file.close();
  } else if (compress == zstd) {
    // not implement yet
  } else {
    std::cout << "type not support yet" << std::endl;
  }

  this->xl = this->boxXH - this->boxXL;
  this->yl = this->boxYH - this->boxYL;
  this->zl = this->boxZH - this->boxZL;
#if DebugInput
  std::cout << "reading complete" << std::endl;
#endif
}

// For the 2d data, by default xy plane is used
void Frame::describe() {
  std::string base = "<FrameData>";
  base += "\n\ttotal particle number= " + std::to_string(particleN) + "\n\tframe time = " + std::to_string(timestep) + "\n<Box size>\n" + "\tx: (" + std::to_string(boxXL) + ", " + std::to_string(boxXH) + ")\n" + "\ty: (" + std::to_string(boxYL) + ", " + std::to_string(boxYH) + ")\n" + "\tz: (" + std::to_string(boxZL) + ", " + std::to_string(boxZH) + ")\n";

  std::string data = "few lines from data:\n";

  if (is2D()) {
    std::cout << base << data
              << std::left << std::setw(12) << "id"
              << std::left << std::setw(12) << "x"
              << std::left << std::setw(12) << "y"
              << "\n";
    for (size_t i = 0; i < std::min(size_t(5), this->particleN); i++)
      std::cout
          << std::left << std::setw(12) << id[i]
          << std::left << std::setw(12) << x[i]
          << std::left << std::setw(12) << y[i]
          << "\n";
    for (size_t i = 0; i < 3; i++)
      std::cout << std::left << std::setw(12) << "...";
    std::cout << "\n";
    for (size_t i = particleN - 5; i < particleN; i++)
      std::cout
          << std::left << std::setw(12) << id[i]
          << std::left << std::setw(12) << x[i]
          << std::left << std::setw(12) << y[i]
          << "\n";
  } else {
    std::cout << base << data
              << std::left << std::setw(12) << "id"
              << std::left << std::setw(12) << "x"
              << std::left << std::setw(12) << "y"
              << std::left << std::setw(12) << "z"
              << "\n";
    for (size_t i = 0; i < std::min(size_t(5), this->particleN); i++)
      std::cout
          << std::left << std::setw(12) << id[i]
          << std::left << std::setw(12) << x[i]
          << std::left << std::setw(12) << y[i]
          << std::left << std::setw(12) << z[i]
          << "\n";
    for (size_t i = 0; i < 4; i++)
      std::cout << std::left << std::setw(12) << "...";
    std::cout << "\n";
    for (size_t i = particleN - 5; i < particleN; i++)
      std::cout
          << std::left << std::setw(12) << id[i]
          << std::left << std::setw(12) << x[i]
          << std::left << std::setw(12) << y[i]
          << std::left << std::setw(12) << z[i]
          << "\n";
  }

  if (!this->attr_order.empty()) {
    std::stringstream ss(this->attr_order);
    std::string token;
    std::cout << "data attributes are:\n";

    while (ss >> token)
      std::cout << std::left << std::setw(12) << token;
    std::cout << "\n";
    for (size_t i = 0; i < std::min(size_t(5), this->particleN); i++) {
      for (size_t j = 0; j < this->attr_table.cols(); j++) {
        std::cout << std::left << std::setw(12) << this->attr_table.data[j][i];
      }
      std::cout << "\n";
    }

    for (size_t i = 0; i < this->attr_table.cols(); i++)
      std::cout << std::left << std::setw(12) << "...";
    std::cout << "\n";

    for (size_t i = particleN - 5; i < particleN; i++) {
      for (size_t j = 0; j < this->attr_table.cols(); j++) {
        std::cout << std::left << std::setw(12) << this->attr_table.data[j][i];
      }
      std::cout << "\n";
    }
  }
}

bool Frame::isInBox() {
  for (size_t i = 0; i < particleN; i++) {
    if ((x[i] < boxXH) && (x[i] > boxXL) && (y[i] < boxYH) && (y[i] > boxYL) && (z[i] < boxZH) && (z[i] > boxZL))
      continue;
    else
      return false;
  }
  return true;
}

bool Frame::is2D() const {
  if (this->z.size() == 0)
    return true;
  return false;
}

bool Frame::is3D() const {
  if (this->z.size() == 0)
    return false;
  return true;
}

void Frame::select(std::string s) {
  // string should be somethin like:
  // c_ackland eq 2
  // the binary symbol follow latex
  // >= geq GEQ
  // == eq EQ
  // <= leq LEQ
  // > g G
  // < l L
  // != neq  NEQ

  const double eps = 1e-8;

  std::vector<size_t> deleteList;
  std::vector<dtype>* target = nullptr;
  std::vector<size_t>* target_int = nullptr;

  std::stringstream ss(s);
  std::string token;
  ss >> token;

  std::function<bool(double)> func;

  if (this->attr_index.count(token)) {
    target = &this->attr_table.data[attr_index[token]];
  } else if (token == "x") {
    target = &this->x;
  } else if (token == "y") {
    target = &this->y;
  } else if (token == "z") {
    target = &this->z;
  } else if (token == "id") {
    target_int = &this->id;
  } else if (token == "type") {
    target_int = &this->type;
  } else {
    throw std::runtime_error("attribute " + token + " not in data!");
  }

  ss >> token;
  double threshold;
  ss >> threshold;

  // eps version
  // if ((token == "geq") || (token == "GEQ") || (token == ">=")) {
  //   func = [threshold, eps](double a) { return a - threshold >= eps; };
  // } else if ((token == "eq") || (token == "EQ") || (token == "==")) {
  //   func = [threshold, eps](double a) { return std::abs(a - threshold) < eps; };
  // } else if ((token == "leq") || (token == "LEQ") || (token == "<=")) {
  //   func = [threshold, eps](double a) { return a - threshold <= eps; };
  // } else if ((token == "g") || (token == "G") || (token == ">")) {
  //   func = [threshold, eps](double a) { return a - threshold > eps; };
  // } else if ((token == "l") || (token == "L") || (token == "<")) {
  //   func = [threshold, eps](double a) { return a - threshold < eps; };
  // } else if ((token == "neq") || (token == "NEQ") || (token == "!=")) {
  //   func = [threshold, eps](double a) { return std::abs(a - threshold) > eps; };
  // }

  // normal version
  if ((token == "geq") || (token == "GEQ") || (token == ">=")) {
    func = [threshold, eps](double a) { return a >= threshold; };
  } else if ((token == "eq") || (token == "EQ") || (token == "==")) {
    func = [threshold, eps](double a) { return a == threshold; };
  } else if ((token == "leq") || (token == "LEQ") || (token == "<=")) {
    func = [threshold, eps](double a) { return a <= threshold; };
  } else if ((token == "g") || (token == "G") || (token == ">")) {
    func = [threshold, eps](double a) { return a > threshold; };
  } else if ((token == "l") || (token == "L") || (token == "<")) {
    func = [threshold, eps](double a) { return a < threshold; };
  } else if ((token == "neq") || (token == "NEQ") || (token == "!=")) {
    func = [threshold, eps](double a) { return a != threshold; };
  }

  if (target != nullptr) {
    for (size_t i = 0; i < target->size(); i++)
      if (!func(target->operator[](i))) {
        // std::cout << target->operator[](i) << token << threshold << std::endl;
        deleteList.push_back(i);
      }
  } else if (target_int != nullptr) {
    for (size_t i = 0; i < target->size(); i++)
      if (!func(target_int->operator[](i)))
        deleteList.push_back(i);
  }
  std::cout << "condition: " << s << std::endl;
  std::cout << "delete " << deleteList.size() << " element of " << this->particleN << std::endl;

  removeIndicesFromVector<dtype, std::vector<size_t>>(this->x, deleteList);
  removeIndicesFromVector<dtype, std::vector<size_t>>(this->y, deleteList);
  removeIndicesFromVector<dtype, std::vector<size_t>>(this->z, deleteList);
  removeIndicesFromVector<size_t, std::vector<size_t>>(this->id, deleteList);
  removeIndicesFromVector<size_t, std::vector<size_t>>(this->type, deleteList);

  for (auto& v : this->attr_table.data)
    removeIndicesFromVector<dtype, std::vector<size_t>>(v, deleteList);

  this->particleN = this->x.size();
}
