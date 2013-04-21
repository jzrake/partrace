#include <iostream>
#include <cmath>
#include "partrace.hpp"

std::ostream &operator<<(std::ostream &s, ChargedParticle &p)
{
  char line[2048];
  snprintf(line, 2048,
	   "%+16.12e %+16.12e %+16.12e %+16.12e "
	   "%+16.12e %+16.12e %+16.12e %+16.12e "
	   "%+16.12e %+16.12e %+16.12e %+16.12e "
	   "%+16.12e %+16.12e %+16.12e %+16.12e",
	   p.x[0], p.x[1], p.x[2], p.x[3],
	   p.u[0], p.u[1], p.u[2], p.u[3],
	   p.B[0], p.B[1], p.B[2], p.B[3],
	   p.E[0], p.E[1], p.E[2], p.E[3]);
  s << line;
  return s;
}

ChargedParticle::ChargedParticle()
{
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  u[0] = 1.0;
  u[1] = 0.0;
  u[2] = 0.0;
  u[3] = 0.0;
  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;
  B[3] = 0.0;
  E[0] = 0.0;
  E[1] = 0.0;
  E[2] = 0.0;
  E[3] = 0.0;
  e = 1.0;
  m = 1.0;
  I = 0;
}

void ParticleManager::print_particle(ChargedParticle &p)
{
  printf("double x[] = { %e, %e, %e, %e };\n", p.x[0], p.x[1], p.x[2], p.x[3]);
  printf("double u[] = { %e, %e, %e, %e };\n", p.u[0], p.u[1], p.u[2], p.u[3]);
  printf("double B[] = { %e, %e, %e, %e };\n", p.B[0], p.B[1], p.B[2], p.B[3]);
  printf("double E[] = { %e, %e, %e, %e };\n", p.E[0], p.E[1], p.E[2], p.E[3]);
}
double ParticleManager::v_dot_v(const double *v)
{
  return v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
}
double ParticleManager::a_dot_b(const double *a, const double *b)
{
  return a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}
void ParticleManager::a_cross_b(const double *a, const double *b, double *c)
{
  c[1] = (a[2]*b[3] - a[3]*b[2]);
  c[2] = (a[3]*b[1] - a[1]*b[3]);
  c[3] = (a[1]*b[2] - a[2]*b[1]);
}


double ChargedParticle::lorentz_factor() const
/* Return the particle's Lorentz factor */
{
  typedef ParticleManager P;
  return sqrt(1 + P::v_dot_v(u));
}
double ChargedParticle::pitch_angle() const
/* Return the cosine of the pitch angle between the velocity and magnetic
   field */
{
  typedef ParticleManager P;
  return P::a_dot_b(u, B) / sqrt(P::v_dot_v(u) * P::v_dot_v(B));
}
double ChargedParticle::larmor_frequency() const
{
  typedef ParticleManager P;
  double B0 = sqrt(P::v_dot_v(B));
  return (e * B0) / (u[0] * m);
}



UniformElectromagneticField::UniformElectromagneticField()
{
  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;
  B[3] = 0.0;
  E[0] = 0.0;
  E[1] = 0.0;
  E[2] = 0.0;
  E[3] = 0.0;
}
void UniformElectromagneticField::
sample_field(ChargedParticle &p) const
{
  std::memcpy(p.B, B, 4 * sizeof(double));
  std::memcpy(p.E, E, 4 * sizeof(double));
}
void UniformElectromagneticField::
set_field(const double B_[4], const double E_[4])
{
  if (B_) std::memcpy(B, B_, 4 * sizeof(double));
  if (E_) std::memcpy(E, E_, 4 * sizeof(double));
}

AlfvenWaveElectromagneticField::AlfvenWaveElectromagneticField()
  : alfven_angular_frequency(1.0)
{

}

void AlfvenWaveElectromagneticField::
sample_field(ChargedParticle &p) const
{
  double B0 = 1.0; // Background field strength
  double B1 = 1e-1; // Field fluctuation
  double va = 0.01; // Alfven speed (in units of c)
  double w = alfven_angular_frequency;
  double k = w / va; // Wavenumber

  double t = p.x[0];
  double x = p.x[1];

  p.B[1] = B0;
  p.B[2] = B1 * cos(k*x - w*t);
  p.B[3] = 0.0;

  p.E[1] = 0.0;
  p.E[2] = 0.0;
  p.E[3] = -va * B1;
}
void AlfvenWaveElectromagneticField::
set_alfven_angular_frequency(double w)
{
  alfven_angular_frequency = 1.0;
}


int compare_movers()
{
  BorisParticleManager boris;
  HymanParticleManager hyman;
  ChargedParticle p1, p2;

  p1.B[0] = 1.0;
  p1.E[1] = 0.5;

  p2.B[0] = 1.0;
  p2.E[1] = 0.5;

  boris.move_particle(p1, 0.01);
  hyman.move_particle(p2, 0.01);

  boris.print_particle(p1);
  hyman.print_particle(p2);

  return 0;
}

int evolve_single()
{
  HymanParticleManager hyman;
  ChargedParticle p1;
  AlfvenWaveElectromagneticField field;

  /*
  UniformElectromagneticField field;
  double B0[4] = { 0, 0, 0, 1 };
  field.set_field(B0, NULL);
  */

  p1.u[2] = 1.0;
  p1.u[0] = p1.lorentz_factor();

  field.sample_field(p1);

  int steps_per_orbit = 100;
  double w = p1.larmor_frequency();
  double f = w / (2 * M_PI);
  double dt = 1.0 / (steps_per_orbit * f);
  double tmax = 1000.0 / f;

  field.set_alfven_angular_frequency(1e-1 * w);

  while (p1.x[0] < tmax) {
    field.sample_field(p1);
    hyman.move_particle(p1, dt);
    std::cout << p1 << std::endl;
  }
  return 0;
}

int main()
{
  evolve_single();
  return 0;
}
