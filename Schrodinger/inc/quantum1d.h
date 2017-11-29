#pragma once

#include "MetComp.h"
//#include "complex.h"
#include <complex>
#include <fftw3.h>

typedef std::complex<real> complex;

/*    Quantum 1D particle    */

class quantum1d {

  real _t;                   // Time
  size_t _N;                 // Space discretisation
  real _dx;                  // Space step
  std::vector<real> _x;      // Position
  std::vector<real> _V;      // Potential
  std::vector<complex> _psi; // Amplitude

  std::function<complex(real real)> _psi_th;

  // static const real _hbar = 6.62607004081e-34;
  // static const real _m = 9.1093835611e-31;
  static const real _hbar = 1;
  static const real _m = 1;

  // Apply normalized hamiltonian (H/ ihbar) to given state
  std::vector<complex> _Hpsi(std::vector<complex> &);

public:
  // Constructors
  template <typename... Ts>
  quantum1d(real a, real b, size_t N, complex (*psi_0_)(real, Ts...),
            Ts... args)
      : _t(0.0), _N(N), _x(std::vector<real>(N)), _V(std::vector<real>(N)),
        _psi(std::vector<complex>(N)) {
    if (a > b)
      std::swap(a, b);

    // Auxiliary variables
    size_t i;
    real x, s;

    // Create space grid
    _dx = (b - a) / (N - 1);

    // Initial wavefunction
    std::function<complex(real)> psi_0 = [psi_0_, args...](real x) {
      return psi_0_(x, args...);
    };

    // Fill vectors
    _x[0] = a;
    for (i = 1; i < _N - 1; i++) {
      x = a + i * _dx;
      _x[i] = x;
      _psi[i] = psi_0(x);
    }
    _x[_N - 1] = b;

    // Normalize
    s = std::sqrt(get_norm());
    for (i = 0; i < _N; i++)
      _psi[i] = _psi[i] / s;
  }

  // Set potential
  template <typename... Ts>
  void set_potential(real (*V_)(real, Ts...), Ts... args) {

    // Auxiliary variables
    size_t i;
    real x;

    // Potential
    std::function<real(real)> V = [V_, args...](real x) {
      return V_(x, args...);
    };

    // Fill vector
    for (i = 0; i < _N; i++)
      _V[i] = V(_x[i]);
  }

  // Update
  void update_euler(real);

  void update_rk4(real);

  void update_btcs(real);

  void update_CN(real);

  void update_euler_FFT(void);

  // Get norm
  real get_norm(void);

  // Plot real part of wavefunction
  void plot_real(void);

  // Plot imaginary part of wavefunction
  void plot_imag(void);

  // Plot probability
  void plot_prob(void);

  // Debug function
  void debug(void);
};
