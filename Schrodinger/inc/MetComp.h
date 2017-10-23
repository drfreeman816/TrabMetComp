#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

typedef double real;

/* Support functions

  print:
  Print series of arguments

  println:
  Print series of arguments with a new line at the end

  printf:
  Print formated text

  printfln:
  Print formated text with a new line at the end

  */
namespace supp {

// Print
template <typename... Ts> void printstr(Ts &&... strs) {
  (std::cout << ... << strs);
}

// Println
template <typename... Ts> void printstrln(Ts &&... strs) {
  (std::cout << ... << strs) << '\n';
}

// Printf
void printf(const char *);

template <typename T, typename... Ts>
void printf(const char *format, T value, Ts... values) {
  for (; *format != '\0'; format++) {
    if (*format == '$') {
      std::cout << value;
      printf(format + 1, values...); // recursive call
      return;
    }
    std::cout << *format;
  }
}

// Printfln
void printfln(const char *);

template <typename T, typename... Ts>
void printfln(const char *format, T value, Ts... values) {
  for (; *format != '\0'; format++) {
    if (*format == '$') {
      std::cout << value;
      printfln(format + 1, values...); // recursive call
      return;
    }
    std::cout << *format;
  }
}

} // namespace supp

/* General math functions

  padder:
  Sum of arguments to a power p

  pnorm:
  Returns p-norm using arguments

  */
namespace math {

static const real pi = std::atan2(0, -1);

// p-Adder
template <typename T> T padder(size_t p, T n) { return pow(n, p); }
template <typename T, typename... Ts> T padder(size_t p, T n, Ts... ns) {
  return pow(n, p) + padder(p, ns...);
}

// p-norm
template <typename... Ts> real pnorm(size_t p, Ts... ns) {
  return pow(padder(p, ns...), 1.0 / p);
}

} // namespace math

/* Complex number */

struct complex {

  bool imag;   // comp: this number has imaginary part
  real re, im; // re: real part, im: imaginary part

  // Constructors
  // Default constructor
  complex(void);
  // Both real and imaginary parts given
  complex(real, real);
  // Only phase given
  complex(real);
  // Real modulus and complex phase given
  complex(real, complex);
  // Complex "modulus" and real phase given (complex + rotation)
  complex(complex, real);

  // Destructor
  ~complex();

  // Swap function
  void swap(complex &);

  // Operators overloading
  // Assignment from a real number
  complex &operator=(const real &);
  // Assignment
  complex &operator=(const complex &);
  // Unary minus
  complex operator-(void) const;
  // Addition with a real number
  complex operator+(const real &) const;
  // Addition
  complex operator+(const complex &) const;
  // Subtraction by a real number
  complex operator-(const real &) const;
  // Subtraction
  complex operator-(const complex &) const;
  // Multiplication by a real number
  complex operator*(const real &)const;
  // Multiplication
  complex operator*(const complex &)const;
  // Division by a real number
  complex operator/(const real &) const;
  // Division
  complex operator/(const complex &) const;

  // Conjugate
  complex conj(void) const;

  // Absolute value
  real abs(void) const;

  // Phase
  real phi(void) const;

  // To string
  std::string tostr(void) const;

}; // Complex number

/* Overload real operators */

/* Complex functions */

// Costant imaginary unit
const complex Iu(0, 1);

// Complex exponential
complex exp(complex);

// Complex logarithm
complex log(complex);

// Complex sine
complex sin(complex);

// Complex cosine
complex cos(complex);

// Complex tangent
complex tan(complex);

// Complex functions

/* Real number
Real number structure with value and uncertainty.
*/
// template <typename T, typename... Ts>
struct Real {

  bool unc;     // unc: this number has uncertainty
  real v, u, r; // v: value, u: uncertainty, r: ratio uncertainty/value

  // Constructors
  Real(real, real);
  Real(real);
  // Destructor
  ~Real();

  // Swap function
  void swap(Real &);

  // Operators overloading
  Real &operator=(const Real &);
  Real operator+(const Real &);
  Real operator-(const Real &);
  Real operator*(const Real &);
  Real operator/(const Real &);
  bool operator<(const Real &);
  bool operator>(const Real &);

  // To string
  std::string tostr(void);

}; // Real number

/* 1D domain */

class Line {

  // Number of points in the line
  size_t _n_points;
  // Points on the line
  std::vector<real> _points;

  // Clear
  void _clear(void);

public:
  // Constructors

  // Endpoints given
  Line(real, real);

  // Vector with points given
  Line(std::vector<real> &);

  // List of points given
  template <typename... Ts> Line(real point, Ts... points) {
    try {
      // Number of arguments passed
      const size_t n_points = sizeof...(Ts) + 1;
      // Add first point
      _points.push_back(point);
      _n_points = 1;
      // Add other points
      add_points(points...);
      // Check if line is created
      if (_points.size() < 2)
        throw 0;
      // Check if all points were valid
      if (_points.size() < n_points) {
        std::cout << "Warning: Some points were not created" << std::endl;
      }
    } catch (...) {
      std::cerr << "Error: no line created" << '\n';
    }
  }

  // Adders

  // Add a point
  void add_point(real);

  // Add a list of points
  void add_points(void) { return; } // Final function call
  template <typename... Ts> void add_points(real point, Ts... points) {
    add_point(point);
    add_points(points...);
  }

  // Add a vector of points
  void add_points(std::vector<real> &);

  // Deleters

  // Delete a point
  void delete_point(real);

  // Delete points (list)
  void delete_points(void) { return; } // Final function call
  template <typename... Ts> void delete_points(real point1, Ts... points) {
    delete_point(point1);
    delete_points(points...);
  }

  // Delete points (vevtor)
  void delete_points(std::vector<real> &);

  // Setters

  // Set endpoints
  void set_endpoints(real, real);

  // Set points (list)
  template <typename... Ts> void set_points(real point, Ts... points) {
    _clear();
    try {
      // Number of arguments passed
      const size_t n_points = sizeof...(Ts) + 1;
      // Add first point
      _points.push_back(point);
      _n_points = 1;
      // Add other points
      add_points(points...);
      // Check if line is created
      if (_points.size() < 2)
        throw 0;
      // Check if all points were valid
      if (_points.size() < n_points) {
        std::cout << "Warning: Some points were not created" << std::endl;
      }
    } catch (...) {
      std::cerr << "Error: no line created" << '\n';
    }
  }

  // Set points (vector)
  void set_points(std::vector<real> &);

  // Getters

  real length(void);

  // Debugging

  // Show all points in order for debugging
  void print_points(void);

}; // Line

/* Multidimensional domain */

class Domain {

  // Domain imension
  size_t _dim;
  std::vector<Line> _lines;

public:
  // Constructors

  // One line given
  Domain(Line);

  // Getters

  // Get Dimension
  size_t dim(void);

  // Debugging

  // Print all points in each line
  void print_lines(void);

}; // Multidimensional domain
