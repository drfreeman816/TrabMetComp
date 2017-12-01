#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Error Class
struct Error {
  char *message;
  char *file;
  size_t line;
  Error(char *m, char *f, size_t l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(Error(message, __FILE__, __LINE__));

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
