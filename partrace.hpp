#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <iostream>


class ChargedParticle;
class ElectromagneticField;
class ParticleManager;


class ElectromagneticField
{
public:
  virtual ~ElectromagneticField() { }
  virtual void sample_field(ChargedParticle &p) const;
} ;

class ChargedParticle
{
public:
  ChargedParticle();
  double x[4], u[4], B[4], E[4];
  double e, m;
  int I;
  double lorentz_factor() const;
} ;

class ParticleManager
{
  friend class ChargedParticle;
public:
  ParticleManager() { }
  virtual ~ParticleManager() { }
  virtual bool move_particle(ChargedParticle &p, double dt) = 0;
  static void print_particle(ChargedParticle &p);
protected:
  static double v_dot_v(const double *v);
  static double a_dot_b(const double *a, const double *b);
  static void a_cross_b(const double *a, const double *b, double *c);
} ;



std::ostream &operator<<(std::ostream &s, ChargedParticle &p);


class UniformElectromagneticField
{
private:
  double B[4];
  double E[4];
public:
  UniformElectromagneticField();
  virtual ~UniformElectromagneticField() { }
  virtual void sample_field(ChargedParticle &p) const;
  void set_field(const double B[4], const double E[4]);
} ;


class HymanParticleManager : public ParticleManager
{
private:
  const double ROOTFIND_TOLER;
  const int NEWTON_ITER_MAX;
  const int BISECT_ITER_MAX;

  double I[4][4], F[4][4], F2[4][4], u[4];
  double b1, b2, e2m;
  bool debug;

public:
  HymanParticleManager();
  void set_debug(bool tf) { debug = tf; }
  bool move_particle(ChargedParticle &p, double dt);

private:
  double evaluate_Eqn32x(double zeta, int mu);
  double evaluate_Eqn32u(double zeta, int mu);
  bool solve_zeta_newton(double &z, double dt);
  bool solve_zeta_bisect(double &z, double dt);
  void print_error(const ChargedParticle &p, double dt);
} ;

class BorisParticleManager : public ParticleManager
{
public:
  bool move_particle(ChargedParticle &p, double dt);
} ;

#endif // PARTICLE_HPP
