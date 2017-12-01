void quantum1d::update_CN_periodic(double dt) {

  size_t i;
  complex a(0, - 0.25 * __dt / (_dx * _dx));
  std::vector<complex> B(_N), u_aux(_N), c_new(_N), d_new(_N), psi_next(_N);

  // Fill B
  for (i = 0; i < _N; i ++) B[i] = complex(1.0, 0.5 * _dt * (1.0 / (_dx * _dx) + _V[i]));

  for (i = 0; i < _N; i ++)  u_aux[i] = conj(a) * (_psi[(i-1+L)%L] + _psi[(i+1+L)%L]) + conj(B[i]) * _psi[i];

  c_new[0] = a / b[0];
  for(i = 1; i < _N; i ++) c_new[i] = a / (B[i] - c_new[i-1] * a);

  d_new[0] = u_aux[0] / B[0];
  for(i = 1; i < _N; i++) d_new[i] = (u_aux[i] - d_new[i-1] * a) / (B[i] - c_new[i-1] * a);

  psi_next[_N-1] = d_new[_N-1];
  for(i = _N-2; i >= 0; i --) psi_next[i] = d_new[i] - c_new[i] * psi_next[i+1];
  
  for(i = 0; i < _N; i++) _psi[i] = psi_next[i];
}

void quantum1d::update_CN_bounded(double dt) {

  size_t i;
  complex a(0, - 0.25 * __dt / (_dx * _dx));
  std::vector<complex> B(_N), u_aux(_N), c_new(_N), d_new(_N), psi_next(_N);

  // Fill B
  for (i = 0; i < _N; i ++) B[i] = complex(1.0, 0.5 * _dt * (1.0 / (_dx * _dx) + _V[i]));

  for (i = 1; i < _N - 1; i ++)  u_aux[i] = conj(a) * (_psi[(i-1+L)%L] + _psi[(i+1+L)%L]) + conj(B[i]) * _psi[i];

  c_new[1] = a / b[1];
  for(i = 2; i < _N - 1; i ++) c_new[i] = a / (B[i] - c_new[i-1] * a);

  d_new[1] = u_aux[1] / B[1];
  for(i = 2; i < _N - 1; i++) d_new[i] = (u_aux[i] - d_new[i-1] * a) / (B[i] - c_new[i-1] * a);

  psi_next[0] = psi_next[_N - 1] = 0;
  psi_next[_N-2] = d_new[_N-2];
  for(i = _N-3; i >= 1; i --) psi_next[i] = d_new[i] - c_new[i] * psi_next[i+1];
  
  for(i = 0; i < _N; i++) _psi[i] = psi_next[i];
}