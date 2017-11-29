#include "MetComp.h"
#include "quantum1d.h"

// Normal distribution1
real normal(real);
// Initial wavefunction
complex psi_0(real);
// Test potential
real V(real);

int main() {
  real a = -10, b = 10;
  size_t N = 1001;
  real dt = 1e-5;

  // Create myparticle
  quantum1d particle(a, b, N, psi_0);

  // Set Potential
  particle.set_potential(V);

  // DEBUG
  particle.debug();

  // Infinite loop
  while (1) {
    std::cerr << "P = " << particle.get_norm() << '\n';
    particle.plot_prob();
    particle.update_rk4(dt);
  }

  return 0;
}

real normal(real x, real mu, real sigma) {
  real N = std::exp(-std::pow(x - mu, 2) / (2 * std::pow(sigma, 2)));
  return (N / (std::sqrt(2 * math::pi) * sigma));
}

complex psi_0(real x) {
  return complex(normal(x, 5, 1), 0) * std::exp(complex(0, math::pi / 4));
}

real V(real x) { return std::pow(x, 2.0); }
