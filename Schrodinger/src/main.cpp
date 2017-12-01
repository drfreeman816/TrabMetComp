#include "MetComp.h"
#include "quantum1d.h"

// Normal distribution1
complex gauss_packet(real, real, real, real);
// Step potential
real porta(real, real, real, real);
// Initial wavefunction
complex psi_0(real);
// Test potential
real V(real);

int main() {
  size_t i;
  real a = -10, b = 10;
  size_t N = 10001;
  real T = 8e-3, dt = 2e-5;
  size_t nt = std::trunc(T / dt);

  std::cerr << nt << '\n';

  // Output files
  std::ofstream norm_dat;

  norm_dat.open("norm.dat");

  // Create myparticle
  quantum1d particle(a, b, N, psi_0);

  // Set Potential
  particle.set_potential(V);

  // DEBUG
  particle.debug();

  // Setup gif
  std::cout << "set terminal gif animate delay 0.1\nset output\'test.gif\'\n";
  std::cout << "set yrange [" << 0 << ':' << 1 << ']' << std::endl;

  // Loop
  // while (1) {
  for (i = 0; i < nt; i++) {
    // Write norm
    norm_dat << nt * dt << '\t' << particle.get_norm() - 1 << '\n';

    particle.plot_prob();
    particle.update_fourier(dt);
    // std::cin.get();
  }

  // Close files
  norm_dat.close();

  return 0;
}

complex gauss_packet(real x, real x_0, real a, real k_0) {
  real A = std::pow(2.0 / (a * math::pi), 0.25);
  return (A * std::exp(complex(-std::pow(x - x_0, 2) / a, k_0 * x)));
}

real porta(real x, real a, real b, real h) {
  if (a > b)
    std::swap(a, b);
  if ((x > a) && (x < b)) {
    return h;
  } else
    return 0;
}

complex psi_0(real x) { return gauss_packet(x, -5, 2, 10); }

real V(real x) { return 1000 * x * x + porta(x, -0.5, 0.5, 10000); }
