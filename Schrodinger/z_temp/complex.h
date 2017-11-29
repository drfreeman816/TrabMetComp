#pragma once

#include "MetComp.h"

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

  // Absolute value squared
  real abs2(void) const;

  // Phase
  real phi(void) const;

  // To string
  std::string tostr(void) const;
  operator std::string() const;

}; // Complex number

/* Overload real operators */

/* Complex functions */

// Costant imaginary unit
const complex Iu(0, 1);

// Square root of a real number
// complex sqrt(real);

// Complex square root
complex sqrt(complex);

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
