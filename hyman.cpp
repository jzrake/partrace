#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "partrace.hpp"

HymanParticleManager::HymanParticleManager() :
  ROOTFIND_TOLER(1e-8),
  NEWTON_ITER_MAX(15),
  BISECT_ITER_MAX(60),
  debug(false)
{
  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      I[mu][nu] = (double)(mu==nu);
    }
  }
}

bool HymanParticleManager::move_particle(ChargedParticle &p, double dt)
{
  const double F_[4][4] = { // Maxwell Tensor F^\mu_\nu
    {  0.0   ,  p.E[1],  p.E[2],  p.E[3] },
    {  p.E[1],  0.0   ,  p.B[3], -p.B[2] },
    {  p.E[2], -p.B[3],  0.0   ,  p.B[1] },
    {  p.E[3],  p.B[2], -p.B[1],  0.0    }
  }; std::memcpy(F, F_, 16*sizeof(double));

  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      F2[mu][nu] = 0.0;
      for (int la=0; la<4; ++la)
	F2[mu][nu] += F[mu][la]*F[la][nu];
    }
  }

  e2m = p.e/p.m;
  b2  = v_dot_v(p.B)-v_dot_v(p.E);
  b1  = sqrt(b2);
  std::memcpy(u, p.u, 4*sizeof(double));

  double zeta = dt*e2m/u[0];
  bool found = false;

  if (!found) {
    found = solve_zeta_newton(zeta,dt);
  }
  if (!found) {
    found = solve_zeta_bisect(zeta,dt);
  }
  if (!found) {
    print_error(p,dt);
    return false;
  }
  else {

    for (int mu=0; mu<4; ++mu)
      {
        p.x[mu] += evaluate_Eqn32x(zeta, mu);
        p.u[mu]  = evaluate_Eqn32u(zeta, mu);
      }
    return true;
  }
}


double HymanParticleManager::evaluate_Eqn32x(double zeta, int mu)
{
  double x_mu=0;
  for (int la=0; la<4; ++la)
    {
      x_mu += u[la] * (I [mu][la] * (b2   * zeta           ) +
                       F [mu][la] * (1    - cos(b1*zeta)   ) +
                       F2[mu][la] * (zeta - sin(b1*zeta)/b1) ) / (e2m*b2);
    }
  return x_mu;
}
double HymanParticleManager::evaluate_Eqn32u(double zeta, int mu)
{
  double u_mu = 0.0;
  for (int la=0; la<4; ++la)
    {
      u_mu += u[la] * (I [mu][la] * (b2                   ) +
                       F [mu][la] * (      sin(b1*zeta)*b1) +
                       F2[mu][la] * (1   - cos(b1*zeta) ) ) / b2;
    }
  return u_mu;
}

bool HymanParticleManager::solve_zeta_newton(double &z, double dt)
{
  double f,g,zeta=z;
  int n_iter = 0;
  do
    {
      f = evaluate_Eqn32x(zeta, 0) - dt;
      g = evaluate_Eqn32u(zeta, 0) / e2m;
      zeta -= f/g;

      if (++n_iter >= NEWTON_ITER_MAX) return false;

    } while (fabs(f) > ROOTFIND_TOLER);
  z = zeta;
  return true;
}
bool HymanParticleManager::solve_zeta_bisect(double &z, double dt)
{
  double f,zeta=z;
  double zetaA = zeta;
  double zetaB = zeta;
  int n_iter = 0;

  while (evaluate_Eqn32x(zetaA, 0) > dt) zetaA *= 0.5;
  while (evaluate_Eqn32x(zetaB, 0) < dt) zetaB *= 2.0;

  zeta = 0.5*(zetaA + zetaB);
  do
    {
      if ((f=evaluate_Eqn32x(zeta, 0) - dt) < 0.0)  zetaA = zeta;
      else                                          zetaB = zeta;
      zeta = 0.5*(zetaA + zetaB);

      if (++n_iter >= BISECT_ITER_MAX) return false;

    } while (fabs(f) > ROOTFIND_TOLER);
  z = zeta;
  return true;
}

void HymanParticleManager::print_error(const ChargedParticle &p, double dt)
{
  printf("## Fatal error in hyman particle mover: rootfinder failed\n");
  printf("dt = %e\n", dt);
  printf("double x[] = { %e, %e, %e, %e };\n", p.x[0], p.x[1], p.x[2], p.x[3]);
  printf("double u[] = { %e, %e, %e, %e };\n", p.u[0], p.u[1], p.u[2], p.u[3]);
  printf("double B[] = { %e, %e, %e, %e };\n", p.B[0], p.B[1], p.B[2], p.B[3]);
  printf("double E[] = { %e, %e, %e, %e };\n", p.E[0], p.E[1], p.E[2], p.E[3]);
}
