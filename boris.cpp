
#include <cmath>
#include "partrace.hpp"

/* ================= Boris Relativistic Particle Mover ================== */
/* ---------------------------------------------------------------------- */
bool BorisParticleManager::move_particle(ChargedParticle &p, double dt)
{
  /* --------------------------------------------------------
   * Relativistic particle mover, attributed to Boris (1970).
   *
   * The algorithm is described in Birdsall & Langdon (1991)
   * Plasma Physics via Computer Simulation, Section 15-4.
   * It is based on a rotation of the (spatial) 4-velocity
   * components around an axis parallel to the magnetic field
   * through an angle
   *
   * theta = -2 arctan( (e B dt) / (2 gamma m c) )
   *
   * The method uses the time-centered 4-velocity, such that
   * p->E and p->B are known on integral time steps t^{n},
   * while p->u is computed for time steps t^{n+1/2}. Note
   * that the positions are known at integral time steps:
   *
   * x^{n+1} = x^{n} + v^{n+1/2} * dt
   *
   * Author: Jonathan Zrake, zrake@nyu.edu
   *
   * --------------------------------------------------------
   */
  double t[4], s[4];
  double u_nph[4], u_nmh[4];
  double u_plus[4], u_minus[4], u_prime[4];
  double u_prime_cross_s[4], u_minus_cross_t[4];
  double gamma_n, gamma_nph;

  double h = (p.e/p.m) * dt;

  for (int d=1; d<=3; ++d)
    u_nmh[d] = p.u[d];

  for (int d=1; d<=3; ++d)                        /* Step 1 */
    u_minus[d] = u_nmh[d] + p.E[d] * 0.5*h;

  gamma_n = sqrt( 1.0 + v_dot_v(u_minus) );       /* Step 2 */

  for (int d=1; d<=3; ++d)                        /* Step 3 */
    t[d] = 0.5*h * p.B[d] / gamma_n;

  for (int d=1; d<=3; ++d)                        /* Step 4 */
    s[d] = 2.0*t[d] / ( 1.0 + v_dot_v(t) );

  a_cross_b(u_minus, t, u_minus_cross_t);

  for (int d=1; d<=3; ++d)                        /* Step 5 */
    u_prime[d] = u_minus[d] + u_minus_cross_t[d];

  a_cross_b(u_prime, s, u_prime_cross_s);

  for (int d=1; d<=3; ++d)                        /* Step 6 */
    u_plus[d] = u_minus[d] + u_prime_cross_s[d];

  for (int d=1; d<=3; ++d)                        /* Step 7 */
    u_nph[d] = u_plus[d] + p.E[d] * 0.5*h;

  gamma_nph = sqrt( 1.0 + v_dot_v(u_nph) );

  p.u[0] = gamma_nph;
  for (int d=1; d<=3; ++d)                        /* Update the velocity */
    p.u[d] = u_nph[d];

  p.x[0] += dt;
  for (int d=1; d<=3; ++d)                        /* Update the position */
    p.x[d] += dt * u_nph[d] / gamma_nph;

  return true;
}

