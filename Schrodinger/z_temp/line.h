#pragma once

#include "MetComp.h"

/* Real function
  f : R x Ts... -> T
 */
template <typename T> class Line {

  // Vector of real number
  typedef std::vector<real> R_vec;

  // Point (real, T)
  typedef std::pair<real, T> Point;
  // Map: vector of points
  typedef std::vector<Point> Map;
  // Iterator over finite elements
  typedef typename Map::iterator Elem_it;

  // Input analytical function
  // typedef std::function<T(real, Ts...)> InFunction;

  // Analytical function stored
  typedef std::function<T(real)> Function;
  // Vector of poiters to analytical functions
  typedef std::vector<Function> Func_vec;

  // Analyitcal part
  /*
  struct Analytical {
    size_t _n_functions;
    Func_vec _functions;
    Args _args;
  };
  */

  // Analytical functions
  Func_vec _functions;
  // Vector with slice points (division of regions between analytical functions)
  R_vec _s_points;
  // Step size in which the domain is diveded
  real _step;
  // Finite elements function evaluations
  Map _elements;
  // Iterator over the finite elements
  Elem_it _ptr;

  // Clear
  void _clear(void) {
    _ptr = nullptr;
    _elements.clear();
  }

  // Lower bound
  Elem_it _lower_bound(real value) {
    return std::lower_bound(
        _elements.begin(), _elements.end(), value,
        [](Point const &p, real x) { return (p.first < x); });
  }

  // Upper bound
  Elem_it _upper_bound(real value) {
    return std::upper_bound(
        _elements.begin(), _elements.end(), value,
        [](real x, Point const &p) { return (x < p.first); });
  }

  // Insert s_point on correct location
  size_t _insert_s_point(real value) {}

public:
  // Constructors

  // Endpoints given
  Line(real value1, real value2) {
    try {
      // Test interval
      if (value1 == value2)
        throw 0;
      // Sort
      if (value1 > value2)
        std::swap(value1, value2);

      // Clear analytical function
      Function zero_function = [](real x) { return 0; };
      _functions.pushback(zero_function);

      // Add points
      _elements.push_back(Point(value1, 0.0));
      _elements.push_back(Point(value2, 0.0));

      _step = value2 - value1;

      // Set element pointer
      _ptr = _elements.begin();

    } catch (...) {
      std::cerr << "Error: endpoints are equal\n";
    }
  }

  // Endpoints and number of points given (equally spaced)
  Line(real value1, real value2, size_t n_points) {
    try {
      // Test interval
      if (value1 == value2)
        throw 0;
      // Sort
      if (value1 > value2)
        std::swap(value1, value2);

      size_t i;

      // Clear analytical function
      Function zero_function = [](real x) { return 0; };
      _functions.pushback(zero_function);

      // Set step
      _step = (value2 - value1) / (n_points - 1);

      // Add first point
      _elements.push_back(Point(value1, 0.0));

      // Add other points
      for (i = 1; i < n_points; i++)
        _elements.push_back(Point(value1 + i * _step, 0.0));

      // Set element pointer
      _ptr = _elements.begin();

    } catch (...) {
      std::cerr << "Error: endpoints are equal\n";
    }
  }

  // Endpoints, number of points and analytical function given (equally spaced)
  template <typename... Ts>
  Line(real value1, real value2, size_t n_points, T (*f)(real, Ts...),
       Ts... args) {
    try {
      // Test interval
      if (value1 == value2)
        throw 0;
      // Sort
      if (value1 > value2)
        std::swap(value1, value2);

      size_t i;
      real x_val;
      T y_val;

      // Set analytical function
      Function function = [f, args...](real x) { return f(x, args...); };
      _functions.push_back(function);

      // Set s_points
      _s_points.push_back(value1);
      _s_points.push_back(value2);

      // Set step
      _step = (value2 - value1) / (n_points - 1);

      // Add points
      for (i = 0; i < n_points; i++) {
        x_val = value1 + i * _step;
        y_val = _functions[0](x_val);
        _elements.push_back(Point(x_val, y_val));
      }

      // Set pointer
      _ptr = _elements.begin();

    } catch (...) {
      std::cerr << "Error: endpoints are equal\n";
    }
  }

  // Iteration operators

  // Deference pointer
  real x(void) { return _ptr->first; }
  T y(void) { return _ptr->second; }

  // Value at index
  real x(size_t index) {
    try {
      if (index >= _elements.size())
        throw 0;

      auto pointer = _elements.begin() + index;

      return pointer->first;

    } catch (...) {
      std::cerr << "Error: Out of bounds. Returning value at pointer." << '\n';
      return _ptr->first;
    }
  }
  T y(size_t index) {
    try {
      if (index >= _elements.size())
        throw 0;

      auto pointer = _elements.begin() + index;

      return pointer->second;

    } catch (...) {
      std::cerr << "Error: Out of bounds. Returning value at pointer." << '\n';
      return _ptr->second;
    }
  }

  // Increment pointer
  Line &operator++() {
    if (_ptr != _elements.end() - 1)
      _ptr++;
    else
      std::cout << "Warning: Already at the end of line." << '\n';
    return *this;
  }

  // Decrement pointer
  Line &operator--() {
    if (_ptr != _elements.begin())
      _ptr--;
    else
      std::cout << "Warning: Already at the beginning of line." << '\n';
    return *this;
  }

  // Return point in index
  Point operator[](size_t index) {
    try {
      if (index >= _elements.size())
        throw 0;

      auto pointer = _elements.begin() + index;

      return *pointer;

    } catch (...) {
      std::cerr << "Erro: Out of bounds" << '\n';
      return *_ptr;
    }
  }

  // Jump to index
  void operator()(size_t index) {
    try {
      if (index >= _elements.size())
        throw 0;

      auto pointer = _elements.begin() + index;
      _ptr = pointer;

    } catch (...) {
      std::cerr << "Erro: Out of bounds" << '\n';
    }
  }

  // Setters

  // Set analytical function everywhere
  template <typename... Ts> void set_function(T (*f)(real, Ts...), Ts... args) {
    // Prepare analytical function
    Function function = [f, args...](real x) { return f(x, args...); };
    _functions.clear();
    _functions.push_back(function);
    // Clear slice points
    _s_points.clear();
    _s_points.push_back(_elements.front().first);
    _s_points.push_back(_elements.back().first);
    // Set values
    for (auto &element : _elements)
      element.second = function(element.first);
  }

  // Set analytical function on intervals (std::function)
  void set_function(real value1, real value2, Function function) {
    // Sort
    if (value1 > value2)
      std::swap(value1, value2);

    // Store number of points
    size_t n_points = _elements.size();

    // Auxiliary variables
    size_t i, j;
    real x_val, x_0;
    T y_val;

    // Set s_points and functions

    // Find slice point >= value1
    R_vec::iterator it_s_points =
        std::lower_bound(_s_points.begin(), _s_points.end(), value1);
    size_t pos = it_s_points - _s_points.begin();
    // Put it there
    it_s_points = ++_s_points.insert(it_s_points, value1);
    // Insert function
    _functions.insert(_functions.begin() + pos, function);
    // Slice previous function
    _functions.insert(_functions.begin() + pos + 1, _functions[pos - 1]);
    // Delete everything until value2
    while (it_s_points !=
           std::upper_bound(_s_points.begin(), _s_points.end(), value2)) {
      it_s_points = _s_points.erase(it_s_points);
      _functions.erase(_functions.begin() + pos + 1);
    }
    // Put it there
    it_s_points = _s_points.insert(it_s_points, value2);

    // Set step
    _step = (_s_points.back() - _s_points.front()) / (n_points - 1);

    // Fill elements
    _elements.clear();
    // First Point
    x_0 = _s_points.front();
    y_val = _functions[0](x_0);
    _elements.push_back(Point(x_0, y_val));
    // Other points
    for (i = 1; i < n_points; i++) {
      x_val = x_0 + i * _step;
      j = std::lower_bound(_s_points.begin(), _s_points.end(), x_val) -
          _s_points.begin();
      y_val = _functions[j - 1](x_val);
      _elements.push_back(Point(x_val, y_val));
    }
  }

  // Set analytical function on interval (pointer to function)
  template <typename... Ts>
  void set_function(real value1, real value2, T (*f)(real, Ts...), Ts... args) {

    // Prepare analytical function
    Function function = [f, args...](real x) { return f(x, args...); };
    // Call main function
    set_function(value1, value2, function);
  }

  // Set endpoints
  void set_endpoints(real value1, real value2) {
    try {
      // Test interval
      if (value1 == value2)
        throw 0;
      // Sort
      if (value1 > value2)
        std::swap(value1, value2);

      // Auxialiry variables
      size_t pos, i, j;
      real x_0, x_val;
      T y_val;

      // Store number of points
      size_t n_points = _elements.size();

      // Set s_points

      // Delete first s_point
      _s_points.erase(_s_points.begin());
      // Find position of value1
      R_vec::iterator it_s_points =
          std::lower_bound(_s_points.begin(), _s_points.end(), value1);
      pos = it_s_points - _s_points.begin();
      // Put it there if its not already
      if (*it_s_points != value1)
        it_s_points = _s_points.insert(it_s_points, value1);
      // Delete everything before that
      _s_points.erase(_s_points.begin(), it_s_points);
      _functions.erase(_functions.begin(), _functions.begin() + pos);

      // Delete last s_point
      _s_points.erase(_s_points.end() - 1);
      // Find position of value2
      it_s_points =
          std::upper_bound(_s_points.begin(), _s_points.end(), value2);
      pos = it_s_points - _s_points.begin();
      // Put it there if its not already
      if (*it_s_points != value2)
        it_s_points = _s_points.insert(it_s_points, value2);
      // Delete everything after that
      _s_points.erase(it_s_points + 1, _s_points.end());
      _functions.erase(_functions.begin() + pos, _functions.end());

      // Set step
      _step = (value2 - value1) / (n_points - 1);

      // Fill elements
      _elements.clear();
      // Fist Point
      x_0 = _s_points.front();
      y_val = _functions[0](x_0);
      _elements.push_back(Point(x_0, y_val));
      // Other points
      for (i = 1; i < n_points; i++) {
        x_val = x_0 + i * _step;
        j = std::lower_bound(_s_points.begin(), _s_points.end(), x_val) -
            _s_points.begin() - 1;
        y_val = _functions[j](x_val);
        //  std::cerr << x_val << '\t' << y_val << '\n';
        _elements.push_back(Point(x_val, y_val));
      }
    } catch (...) {
      std::cerr << "Error: endpoints are equal\n";
    }
  }

  void set_n_points(size_t n_points) {
    try {
      if (n_points < 2)
        throw 0;

      // Auxialiry variables
      size_t i, j;
      real x_0, x_val;
      T y_val;

      // Set step
      _step = (_s_points.back() - _s_points.front()) / (n_points - 1);

      // Fill elements
      _elements.clear();
      // Fist Point
      x_0 = _s_points.front();
      y_val = _functions[0](x_0);
      _elements.push_back(Point(x_0, y_val));
      // Other points
      for (i = 1; i < n_points; i++) {
        x_val = x_0 + i * _step;
        j = std::lower_bound(_s_points.begin(), _s_points.end(), x_val) -
            _s_points.begin() - 1;
        y_val = _functions[j](x_val);
        //  std::cerr << x_val << '\t' << y_val << '\n';
        _elements.push_back(Point(x_val, y_val));
      }

    } catch (...) {
      std::cerr << "Error: at least two points are required\n";
    }
  }

  void set_zero(real value1, real value2) {
    Function zero_function = [](real x) -> T { return 0; };
    set_function(value1, value2, zero_function);
  }

  // Getters

  // Analytical function associated
  T f(real x) {

    // Check boundaries
    if (x <= _s_points.back())
      return _functions.back()(x);
    if (x >= _s_points.front())
      return _functions.front()(x);

    R_vec::iterator it_s_points =
        std::lower_bound(_s_points.begin(), _s_points.end(), x);

    size_t j = it_s_points - _s_points.begin() - 1;

    return _functions[j](x);
  }

  // Get start point
  real start(void) { return _elements.front().first; }
  // Get final point
  real end(void) { return _elements.back().first; }

  // Get length
  real length(void) {
    return (_elements.back().first - _elements.front().first);
  }

  // Set Step
  real step(void) { return _step; }

  // Set y-value
  void set_y(size_t index, T value) { _elements[index].second = value; }

  // Get number of points
  size_t n_points(void) { return _elements.size(); }

  // Get number of analytical functions
  size_t n_func(void) { return _functions.size(); }

  // Debugging

  // Show all points in order for debugging
  void print_points(void) {
    std::cerr << "print_points:" << '\n';
    // typename Map::iterator it;
    for (const auto &element : _elements)
      std::cerr << element.first << "\t\t" << std::string(element.second)
                << '\n';
    std::cerr << "\n";
  }

  // Show all s_points in order for debugging
  void print_s_points(void) {
    std::cerr << "print_s_points:" << '\n';
    // typename Map::iterator it;
    for (const auto &s_point : _s_points)
      std::cerr << s_point << '\n';
    std::cerr << "\n";
  }

  // Debug
  void debug(void) {
    std::cerr << "DEBUG LOG\n\n";

    std::cerr << "number of points = " << n_points() << '\n';
    std::cerr << "number of functions = " << n_func() << '\n';

    print_points();
    print_s_points();

    std::cerr << "length = " << length() << "\n\n";
  }
}; // Line
