// Implementation of Kristensen et al 2006 CJFAS

#include <TMB.hpp>

template<class Type>
struct model{
  // Temporary variables
  Type t0_global;
  Type Linf_global;

  // Model parameters
  Type L0;                // Size at recruitment
  Type k;                 // von Bertalanffy k
  Type M0;                // Natural mortality (constant by size)
  Type beta;              // Fishing mortality parameter
  Type l50f;              // Fishing mortality parameter
  Type Finf;              // Fishing mortality parameter
  Type alpha;             // Survey selection parameter
  Type l50;               // Survey selection parameter
  Type meanLinf;          // E[ L_infinity ]
  Type sdLinf;            // SD[ L_infinity ]
  Type p;                 // 'Effort' parameter
  Type minsel;            // Not used. Set to 0
  vector<Type> R;         // Recruitment vector
  vector<Type> t0;        // Recruitment times
  vector<Type> spawnsd;   // SD of recruitment peaks

  // Tuning parameteres for integration accuracy
  int nsteps;   // romberg
  int promberg; // romberg
  int n_sp;
  int maxage;
  // von Bertalanffy
  Type L(Type t) {
    return Linf_global - (Linf_global-L0) * exp( - k * (t - t0_global) );
  }
  // Natural mortality as function of size and time
  Type M(Type x, Type t)
  {
    return M0;
  }
  // Fishing mortality as function of size
  Type Fx(Type x)
  {
    return 1.0 / (1.0 + exp( -beta * (x - l50f) ) );
  }
  // Fishing mortality as function of time
  Type Ft(Type t)
  {
    return Finf;
  }
  // Fishing mortality as function of size and time
  Type F(Type x, Type t){
    return Fx(x) * Ft(t);
  }
  // Total mortality as function of size and time
  Type z(Type x, Type t)
  {
    return F(x,t) + M(x,t);
  }   
  // Integrand (see paper) f
  Type f(Type t)
  {
    return z(L(t), t);
  }
  struct integrand{
    model<Type>* pmod;
    integrand (model<Type>* pmod_) { pmod=pmod_; }
    Type operator()(const Type &t) { return pmod->f(t); }
  };
  Type romberg(Type a, Type b, int n){
    Type e;
    integrand f(this);
    return CppAD::RombergOne(f, a, b, n, promberg, e);
  }
  // Notation as in paper
  Type Z(Type t, Type t0, Type Linf)
  {
    t0_global=t0;
    Linf_global=Linf;
    return romberg(t0,t,log2(nsteps)+2);
  }
  // Density of individual L_infinity
  Type u(Type x)
  {
    return dnorm(x, meanLinf, sdLinf);
  }
  // Survey selection
  Type s(Type x)
  {
    Type ans = 1 / (1 + exp(-alpha*(x-l50)));
    ans = (ans + minsel) / (Type(1) + minsel);
    return ans;
  }
  
  Type q(Type x, int i)
  {
    return p * s(x,i);
  }
  // Recruitment - sum of gaussian peaks
  Type r(Type t)
  {
    Type res = 0.0;
    for(int i=0; i<t0.size(); i++) res += R[i] * dnorm(t, t0[i], spawnsd[i]);
    return res;
  }
  // See paper
  Type n_integrand(Type x, Type t, Type tinp)
  {
    Type t1 = tinp;
    Type tmp = exp(-k*(t-t1));
    Type dGx = 1 / (1 - tmp);
    dGx = CppAD::CondExpLe(t, t1, Type(0.0), dGx);
    Type Gx = (x - L0*tmp) * dGx;
    Type res = r(t1) * exp(-Z(t,t1,Gx)) * u(Gx) * dGx;
    return res;
  }
  /* For Romberg integration */
  class integrand2{
  public:
    model<Type>* pmod;
    Type x,t;
    integrand2 (model<Type>* pmod_, Type x_, Type t_){pmod=pmod_; x=x_; t=t_;}
    Type operator()(const Type &s){return pmod->n_integrand(x,t,s);}
  };
  
  // Number density
  Type n(Type x, Type t)
  {
    Type a = t-maxage;
    Type b = t;
    int n = (n_sp/2)*int(maxage);
    n=log2(n)+2;
    Type e;
    integrand2 f(this,x,t);
    return CppAD::RombergOne(f,a,b,n,promberg,e);
  }
  // Number density times survey selection (To be linked with observed counts)
  Type ns(Type x, Type t){
    return n(x,t) * s(x);
  }

};


template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(L0);
  PARAMETER(meanLinf);
  PARAMETER(sdLinf);
  PARAMETER(k);
  PARAMETER(beta);
  PARAMETER(l50f);
  PARAMETER(M0);
  PARAMETER(l50);
  PARAMETER(alpha);
  PARAMETER(minsel);
  PARAMETER(logp);
  PARAMETER(Finf);
  PARAMETER_VECTOR(logR);
  PARAMETER_VECTOR(t0);
  PARAMETER(dt0);  // Birth date
  t0 = t0 + dt0;
  PARAMETER_VECTOR(spawnsd);
  DATA_VECTOR(times); /* survey times */
  DATA_VECTOR(sizes); /* Observed sizes  */
  DATA_MATRIX(N);  /* Observation matrix - size fastest running: N[size,time] */
  model<Type> mod;
  mod.L0=L0;
  mod.k=k;
  mod.M0=M0;
  mod.beta=beta;
  mod.l50f=l50f;
  mod.Finf=Finf;
  mod.alpha=alpha;
  mod.l50=l50;
  mod.meanLinf=meanLinf;
  mod.sdLinf=sdLinf;
  mod.p=exp(logp);
  // Vector
  mod.R=exp(logR);
  mod.t0=t0;
  mod.spawnsd=spawnsd;
  mod.nsteps=10;   // romberg
  mod.promberg=3; // romberg
  mod.n_sp=10;
  mod.maxage=8;
  mod.minsel=minsel;
  matrix<Type> mat(N);
  for(int i=0; i<times.size(); i++)
    for(int j=0; j<sizes.size(); j++)
      mat(j,i)=mod.ns(sizes[j],times[i]);

  REPORT(mat);
  REPORT(N);

  // Least squares for simplicity (Kristensen et al used negative binomial)
  matrix<Type> diff = mat - N;
  return (diff.array() * diff.array()).sum();
}
