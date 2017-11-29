#include "../inc/quantum1d.h"

/*    Quantum 1D particle    */

// Apply normalized hamiltonian (H/ ihbar) to given state
std::vector<complex> quantum1d::_Hpsi(std::vector<complex> &psi) {
  size_t i;
  std::vector<complex> Hpsi(_N);
  complex z1 = complex(0.0, _hbar / (2 * _m * _dx * _dx)), z2;

  for (i = 1; i < _N - 1; i++) {
    z2 = complex(0.0, -_V[i] / _hbar);
    Hpsi[i] = z1 * (psi[i - 1] - 2 * psi[i] + psi[i + 1]) - z2 * psi[i];
  }

  return Hpsi;
}

// Update
void quantum1d::update_euler(real dt) {
  size_t i;
  real k = dt / (2.0 * _dx * _dx);
  complex z, p0, p1, p2;

  p0 = _psi[0];
  p1 = _psi[1];
  for (i = 1; i < _N - 1; i++) {
    z = complex(1.0, -2.0 * k * _hbar / _m - dt * _V[i] / _hbar);

    p2 = _psi[i + 1];

    _psi[i] = z * p1 + complex(0.0, k * _hbar / _m) * (p0 + p2);

    p0 = p1;
    p1 = p2;
  }

  _t += dt;
}

void quantum1d::update_rk4(real dt) {
  size_t i;
  real dt2 = dt / 2;
  std::vector<complex> psi_k1(_N), psi_k2(_N), psi_k3(_N), psi_k4(_N),
      psi_temp(_N);

  // K1
  psi_k1 = _Hpsi(_psi);
  // K2
  for (i = 1; i < _N - 1; i++)
    psi_temp[i] = _psi[i] + dt2 * psi_k1[i];
  psi_k2 = _Hpsi(psi_temp);
  // K3
  for (i = 1; i < _N - 1; i++)
    psi_temp[i] = _psi[i] + dt2 * psi_k2[i];
  psi_k3 = _Hpsi(psi_temp);
  // K4
  for (i = 1; i < _N - 1; i++)
    psi_temp[i] = _psi[i] + dt * psi_k3[i];
  psi_k4 = _Hpsi(psi_temp);

  // Apply
  for (i = 1; i < _N - 1; i++)
    _psi[i] += dt * (psi_k1[i] + 2 * psi_k2[i] + 2 * psi_k3[i] + psi_k4[i]) / 6;

  _t += dt;
}

// Crank-Nicolson method
void quantum1d::update_CN(double dt) {

  size_t i, n = _N - 2;

  double k = dt / (2 * _dx * _dx);
  complex beta(0, -k / 2);
  std::vector<real> B(n);
  std::vector<complex> alpha(n), chi(n), gamma(n), phi(n);

  // Fill B
  for (i = 0; i < n; i++) {
    B[i] = k + dt * _V[i + 1] / 2;
  }

  // Fill alpha
  for (i = 0; i < n; i++)
    alpha[i] = complex(1.0, B[i]);

  // Fill chi
  for (i = 0; i < n; i++) {
    chi[i] = alpha[i] * _psi[i + 1] - beta * (_psi[i] + _psi[i + 2]);
  }

  // Calculate gamma
  gamma[0] = beta / alpha[0];
  for (i = 1; i < n; i++)
    gamma[i] = beta / (alpha[i] - beta * gamma[i - 1]);

  // Calculate phi
  phi[0] = (chi[0] - beta * _psi[0]) / alpha[0];
  for (i = 1; i < n - 1; i++)
    phi[i] = (chi[i] - beta * phi[i - 1]) / (alpha[i] - beta * gamma[i - 1]);
  phi[n - 1] = (chi[n - 1] - beta * (_psi[_N - 1] + phi[n - 2])) /
               (alpha[n - 1] - beta * gamma[n - 2]);

  // Update psi
  _psi[_N - 2] = phi[n];
  for (i = _N - 3; i > 0; i--) {
    _psi[i] = phi[i - 1] - gamma[i - 1] * _psi[i + 1];
  }
}

void quantum1d::update_euler_FFT(void) {

  size_t i;
  real k;

  // FFTW3 arrays
  fftw_complex *psi, *phi;
  // FFTW3 plans
  fftw_plan ft, ift;

  // Allocate arrays
  psi = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * _N);
  phi = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * _N);

  // Specify FFTW plans
  ft = fftw_plan_dft_1d(_N, psi, phi, FFTW_FORWARD, FFTW_MEASURE);
  ift = fftw_plan_dft_1d(_N, phi, psi, FFTW_BACKWARD, FFTW_MEASURE);

  // Fill psi
  for (i = 0; i < _N; i++)
    psi[i] = _psi[i];

  // Execute FFTW plan
  fftw_execute(ft);

  // Apply diagonal operator
  for (i = 0; i < _N; i++) {
    k = i - _N / 2;
  }

  // Destroy FFTW plans
  fftw_destroy_plan(ft);
  fftw_destroy_plan(ift);
  // Free memory allocated
  fftw_free(psi);
  fftw_free(phi);
}

// Get norm
real quantum1d::get_norm(void) {
  real S = 0;

  for (const auto &A : _psi)
    S += std::norm(A);

  return _dx * S * (_N - 1) / _N;
}

// Plot real part of wavefunction
void quantum1d::plot_real(void) {
  // Setup GNUPLOT
  std::cout << "set key off" << std::endl;
  std::cout << "set xrange [" << _x.front() << ':' << _x.back() << ']'
            << std::endl;

  // Call interactive terminal
  std::cout << "plot \"-\" w p pt 7 ps 0.5" << std::endl;

  for (size_t i = 0; i < _N; i++)
    std::cout << _x[i] << "\t" << _psi[i].real() << '\n';
  std::cout << 'e' << std::endl;
}

// Plot imaginary part of wavefunction
void quantum1d::plot_imag(void) {
  // Setup GNUPLOT
  std::cout << "set key off" << std::endl;
  std::cout << "set xrange [" << _x.front() << ':' << _x.back() << ']'
            << std::endl;

  // Call interactive terminal
  std::cout << "plot \"-\" w p pt 7 ps 0.5" << std::endl;

  for (size_t i = 0; i < _N; i++)
    std::cout << _x[i] << "\t" << _psi[i].imag() << '\n';
  std::cout << 'e' << std::endl;
}

// Plot probability
void quantum1d::plot_prob(void) {
  // Setup GNUPLOT
  std::cout << "set key off" << std::endl;
  std::cout << "set xrange [" << _x.front() << ':' << _x.back() << ']'
            << std::endl;

  // Call interactive terminal
  std::cout << "plot \"-\" w p pt 7 ps 0.5" << std::endl;

  for (size_t i = 0; i < _N; i++)
    std::cout << _x[i] << "\t" << std::norm(_psi[i]) << '\n';
  std::cout << 'e' << std::endl;
}

// Debug function
void quantum1d::debug(void) {
  std::cerr << "1D quantum particle debug\n\n";

  std::cerr << "Domain = [" << _x.front() << ", " << _x.back() << "]\n";
  std::cerr << "N = " << _N << " -> dx = " << _dx << "\n\n";

  std::cerr << "_x.size() = " << _x.size() << '\n';
  std::cerr << "_V.size() = " << _V.size() << '\n';
  std::cerr << "_psi.size() = " << _psi.size() << "\n\n";

  std::cerr << "P_total = " << get_norm() << "\n\n";

  std::cerr << "t = " << _t << "\n\n";

  std::cerr << "x\t\tV(t,x)\t\tpsi(t,x)" << '\n';
  for (size_t i = 0; i < _N; i++)
    std::cerr << _x[i] << "\t\t" << _V[i] << "\t\t" << _psi[i] << '\n';
}
