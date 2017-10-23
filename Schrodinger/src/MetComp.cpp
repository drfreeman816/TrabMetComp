#include "../inc/MetComp.h"

/* Support functions */

namespace supp {
// Printf
void printf(const char *format) { std::cout << format; }

// Printfln
void printfln(const char *format) { std::cout << format << std::endl; }
} // namespace supp

/* Math */

// namespace math

/* Complex number */

// Constructors

// Default constructor
complex::complex(void) : imag(false), re(0), im(0) {}

// Both real and imaginary parts given
complex::complex(real re0, real im0) : imag(true), re(re0), im(im0) {
  if (im0 == 0)
    imag = false;
}

// Only phase given
complex::complex(real phi0)
    : imag(true), re(std::cos(phi0)), im(std::sin(phi0)) {
  if (im == 0)
    imag = false;
}

// Real modulus and complex phase given
complex::complex(real r0, complex phi0) : imag(true) {
  *this = exp(phi0) * r0;
  if (im == 0)
    imag = false;
}

// Complex "modulus" and real phase given (complex + rotation)
complex::complex(complex r0, real phi0) : imag(true) {
  *this = r0 * exp(Iu * phi0);
  if (im == 0)
    imag = false;
}

// Destructor
complex::~complex() {}

// Swap function
void complex::swap(complex &other) {
  std::swap(imag, other.imag);
  std::swap(re, other.re);
  std::swap(im, other.im);
}

// Operators overloading

// Assignment from a real number
complex &complex::operator=(const real &other) {
  complex temp(other, 0);
  this->swap(temp); // this <=> temp
  return *this;
}

// Assignment
complex &complex::operator=(const complex &other) {
  complex temp(other); // temp = other
  this->swap(temp);    // this <=> temp
  return *this;
}

// Unary minus
complex complex::operator-(void) const {
  complex result(-re, -im);
  return result;
}

// Addition with a real number
complex complex::operator+(const real &other) const {
  complex result(re + other, im);
  return result;
}

// Addition
complex complex::operator+(const complex &other) const {
  complex result(re + other.re, im + other.im);
  return result;
}

// Subtraction by a real number
complex complex::operator-(const real &other) const {
  complex result(re - other, im);
  return result;
}

// Subtraction
complex complex::operator-(const complex &other) const {
  complex result(re - other.re, im - other.im);
  return result;
}

// Multiplication by a real number
complex complex::operator*(const real &other) const {
  complex result(re * other, im * other);
  return result;
}

// Multiplication
complex complex::operator*(const complex &other) const {
  complex result(re * other.re - im * other.im, re * other.im + im * other.re);
  return result;
}

// Division by a real number
complex complex::operator/(const real &other) const {
  complex result(re / other, im / other);
  return result;
}

// Division
complex complex::operator/(const complex &other) const {
  complex result = *this * other.conj();
  result = result / math::padder(2, other.re, other.im);
  return result;
}

// Conjugate
complex complex::conj(void) const {
  complex result(re, -im);
  return result;
}

// Absolute value
real complex::abs(void) const {
  if (!imag)
    return std::abs(re);
  else
    return math::pnorm(2, re, im);
}

// Phase
real complex::phi(void) const {
  if (!imag) {
    if (re > 0)
      return 0;
    else
      return math::pi;
  } else
    return std::atan(re / im);
}

// To string
std::string complex::tostr(void) const {
  if (!imag)
    return std::to_string(re);
  else if (im > 0)
    return (std::to_string(re) + " + " + std::to_string(im) + "i");
  else
    return (std::to_string(re) + " - " + std::to_string(std::abs(im)) + "i");
}

// Complex number

/* Overload real operators */

/* Complex functions */

// Complex exponential
complex exp(complex val) {
  complex result(std::exp(-val.im) * std::cos(val.re),
                 std::exp(-val.im) * std::sin(val.re));
  return result;
}

// Complex logarithm
complex log(complex val) {
  complex result(std::log(val.abs()), val.phi());
  return result;
}

// Complex sine
complex sin(complex val) {
  complex result;
  result = (exp(Iu * val) - exp(-Iu * val)) / 2;
  return result;
}

// Complex cosine
complex cos(complex val) {
  complex result;
  result = (exp(Iu * val) + exp(-Iu * val)) / 2;
  return result;
}

// Complex tangent
complex tan(complex val) {
  complex result;
  result = sin(val) / cos(val);
  return result;
}

// Complex functions

/* Real number */

// Constructors
Real::Real(real v0, real r0)
    : unc(true), v(v0), u(std::abs(r0 * v0)), r(std::abs(r0)) {
  try {
    if (v0 == 0) {
      unc = false;
      u = 0;
      r = 0;
      throw "Uncertainty undefined for value zero.";
    } else if (r0 == 0) {
      unc = false;
      u = 0;
    }
  } catch (const char *msg) {
    std::cerr << msg << std::endl;
  }
}

Real::Real(real v0) : Real(v0, 0) { unc = false; }

// Destructor
Real::~Real() {}

// Swap function
void Real::swap(Real &other) {
  std::swap(unc, other.unc);
  std::swap(v, other.v);
  std::swap(r, other.r);
  std::swap(u, other.u);
}

// Operators overloading

// Assignment
Real &Real::operator=(const Real &other) {
  Real temp(other);
  this->swap(temp);
  return *this;
}

// z = x + y
Real Real::operator+(const Real &y) {
  Real z(v + y.v);
  if (unc || y.unc) {
    z.unc = true;
    z.u = math::pnorm(2, u, y.u);
    z.r = z.u / z.v;
  }
  return z;
}
// z = x - y
Real Real::operator-(const Real &y) {
  Real z(v - y.v);
  if (unc || y.unc) {
    z.unc = true;
    z.u = math::pnorm(2, u, y.u);
    z.r = z.u / z.v;
  }
  return z;
}
// z = x * y
Real Real::operator*(const Real &y) {
  Real z(v * y.v);
  if (unc || y.unc) {
    z.unc = true;
    z.u = math::pnorm(2, v * u, y.v * y.u);
    z.r = z.u / z.v;
  }
  return z;
}
// z = x / y
Real Real::operator/(const Real &y) {
  try {
    if (y.v == 0)
      throw "Cannot divide by zero! What the fu-";
    Real z(v / y.v);
    if (unc || y.unc) {
      z.unc = true;
      z.u = math::pnorm(2, u / y.v, v * y.u / (y.v * y.v));
      z.r = z.u / z.v;
    }
    return z;
  } catch (const char *msg) {
    std::cerr << msg << '\n';
    return *this;
  }
}
// x < y
bool Real::operator<(const Real &y) {
  try {
    real d = u + y.u / (std::abs(v - y.v));
    if (d < 1)
      return (v < y.v);
    else
      throw d;
  } catch (real d) {
    std::cerr << "Uncertainty overlap: d = " << d << '\n';
    return (v < y.v);
  }
}
// x > y
bool Real::operator>(const Real &y) {
  try {
    real d = u + y.u / (std::abs(v - y.v));
    if (d < 1)
      return (v > y.v);
    else
      throw d;
  } catch (real d) {
    std::cerr << "Uncertainty overlap: d = " << d << '\n';
    return (v > y.v);
  }
}

// To string
std::string Real::tostr(void) {
  if (unc)
    return (std::to_string(v) + " +/- " + std::to_string(u) +
            " (r=" + std::to_string(r) + ")");
  else
    return std::to_string(v);
}

/* Line */

// Constructors

// Endpoints given
Line::Line(real point1, real point2) {
  try {
    if (point1 == point2)
      throw 0;
    _points.push_back(point1);
    _n_points = 1;
    add_point(point2);
  } catch (...) {
    std::cerr << "Error: endpoints are equal" << '\n';
  }
}

// Vector with points given
Line::Line(std::vector<real> &points) {
  try {
    // Number of arguments passed
    size_t n_points = points.size();
    // Add first point
    _points.push_back(points.back());
    _n_points = 1;
    points.pop_back();
    // Add other points
    add_points(points);
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

// Clear
void Line::_clear(void) {
  _points.clear();
  _n_points = 0;
}

// Adders

// Add a point
void Line::add_point(real point) {
  // Iterate on points and find pointer to first element larger than point
  std::vector<real>::iterator it_points =
      std::upper_bound(_points.begin(), _points.end(), point);
  // Insert point in the correct place if its not already there
  if (*(it_points - 1) != point) {
    _points.insert(it_points, point);
    _n_points++;
  }
}

// Add a vector of points
void Line::add_points(std::vector<real> &points) {
  // Iterate on vector and add each point
  std::vector<real>::iterator it_points;
  for (it_points = points.begin(); it_points != points.end(); it_points++)
    add_point(*it_points);
}

// Deleters

// Delete a point
void Line::delete_point(real point) {
  try {
    // Are there enough points
    if (_points.size() < 3)
      throw point;
    // Look for given point
    std::vector<real>::iterator it_points =
        std::find(_points.begin(), _points.end(), point);
    // Delete it if this point exists
    if (*it_points == point)
      _points.erase(it_points);
    else
      std::cerr << "Warning: Point " << point << " not defined in this line."
                << '\n';
  } catch (real err_point) {
    std::cerr << "Error: This line only have two points" << '\n';
    std::cerr << "Point " << err_point << " not deleted." << '\n';
    return;
  }
}

// Delete points (list)

// Delete points (vevtor)
void Line::delete_points(std::vector<real> &points) {
  // Iterate over vector and delete each point
  std::vector<real>::iterator it_points = points.begin();
  for (; it_points != points.end(); it_points++)
    delete_point(*it_points);
}

// Setters

// Set endpoints
void Line::set_endpoints(real point1, real point2) {
  std::vector<real>::iterator it_points;
  // Check order of parameters
  if (point2 < point1)
    std::swap(point1, point2);
  // Find lower bound, put point1 there and erase everything before that
  it_points = std::lower_bound(_points.begin(), _points.end(), point1);
  if (*it_points != point1)
    it_points = _points.insert(it_points, point1);
  _points.erase(_points.begin(), it_points);
  // Find upper bound, put point2 there and erase everything after that
  it_points = std::upper_bound(_points.begin(), _points.end(), point2);
  if (*it_points != point2)
    it_points = _points.insert(it_points, point2);
  _points.erase(++it_points, _points.end());
}

// Set points (list)

// Set points (vector)
void Line::set_points(std::vector<real> &points) {
  _clear();
  try {
    // Number of arguments passed
    size_t n_points = points.size();
    // Add first point
    _points.push_back(points.back());
    _n_points = 1;
    points.pop_back();
    // Add other points
    add_points(points);
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

// Getters

real Line::length(void) { return (_points.back() - _points.front()); }

// Debugging

// Show all points in order for debugging
void Line::print_points(void) {
  std::cerr << "\nprint_points" << '\n';
  std::vector<real>::iterator it_points;
  for (it_points = _points.begin(); it_points != _points.end(); it_points++)
    std::cerr << "*it_points = " << *it_points << '\n';
  std::cerr << "\n";
}

// Line

/* Multidimensional domain */

// Constructors

// One line given
Domain::Domain(Line line) : _dim(1) { _lines.push_back(line); }

// Getters

// Get Dimension
size_t Domain::dim(void) { return _dim; }

// Debugging

// Print all points in each line
void Domain::print_lines(void) {
  std::cerr << "\nprint_lines" << '\n';
  std::vector<Line>::iterator it_lines;
  for (it_lines = _lines.begin(); it_lines != _lines.end(); it_lines++)
    it_lines->print_points();
  std::cerr << "\n";
}
