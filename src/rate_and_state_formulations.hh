/* -------------------------------------------------------------------------- */
struct StateEvolution {
  StateEvolution(){};
  virtual ~StateEvolution(){};
  virtual inline Real operator()(Real rate, Real state, Real D, Real delta_t){
    return (1-rate*state/D)*delta_t;
  }
  virtual inline Real getSteadyState(Real rate, Real D) {
    return D/rate;
  }
};

/* -------------------------------------------------------------------------- */
struct SlipStateEvolution : public StateEvolution {
  SlipStateEvolution(){};
  virtual ~SlipStateEvolution(){};
  virtual inline Real operator()(Real rate, Real state, Real D, Real delta_t){
    return -rate*state/D*log(rate*state/D)*delta_t;
  }
  virtual inline Real getSteadyState(Real rate, Real D) {
    return D/rate;
  }
};

/* -------------------------------------------------------------------------- */
struct RegularizedStateEvolution : public StateEvolution {
  RegularizedStateEvolution(Real v0){this->v0=v0;};
  virtual ~RegularizedStateEvolution(){};
  inline Real operator()(Real rate, Real state, Real D, Real delta_t){
    return (1-sqrt(rate*rate+v0*v0)*state/D)*delta_t;
  }
  virtual inline Real getSteadyState(Real rate, Real D) {
    return D/sqrt(rate*rate+v0*v0);
  }

private:
  Real v0;
};

/* -------------------------------------------------------------------------- */
struct RandSFormulation {
  RandSFormulation(){};
  virtual ~RandSFormulation(){};
  virtual inline Real operator()(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return f0 + a*log(1+rate/rate_st)+b*f0*log(1+state/state_st);
  }
  virtual inline Real getTangent(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return a/(rate_st+rate);
  }
  virtual inline Real getSteadyTangent(Real rate, Real state, Real a, Real b, Real D,
				       Real f0, Real rate_st, Real state_st){
    return a/(rate_st+rate) - b*f0*D/(rate*rate*state_st+rate*D);
  }
  virtual inline Real getStableState(Real interface_traction, Real sigma_0, Real rate,
				     Real a, Real b, Real D, Real f0, Real rate_st, Real state_st) {
    return 0;
  }
};

/* -------------------------------------------------------------------------- */
struct WeakeningRandSFormulation : public RandSFormulation {
  WeakeningRandSFormulation(){};
  virtual ~WeakeningRandSFormulation(){};
  virtual inline Real operator()(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return f0 + a*log(rate/rate_st)+b*f0*log(state/state_st);
  }
    
  virtual inline Real getTangent(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return a/rate;
  }
  virtual inline Real getSteadyTangent(Real rate, Real state, Real a, Real b, Real D,
				       Real f0, Real rate_st, Real state_st){
    return (a-b*f0)/rate - b*f0/rate;
  }
  virtual inline Real getStableState(Real interface_traction, Real sigma_0, Real rate,
				     Real a, Real b, Real D, Real f0, Real rate_st, Real state_st) {
    return 0;
  }
};

/* -------------------------------------------------------------------------- */
struct RegularizedRandSFormulation : public RandSFormulation {
  RegularizedRandSFormulation(Real v0, Real theta, Real xi){
    this->v0=v0; this->theta=theta; this->xi=xi;}
  virtual ~RegularizedRandSFormulation(){};
  virtual inline Real operator()(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return (1+b*log(1+state/state_st))*(theta/sqrt(1+v0*v0/(rate*rate))+xi*log(1+rate/rate_st));
  }
  virtual inline Real getTangent(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return (1+b*log(1+state/state_st))*(theta*v0*v0/(pow(rate*rate+v0*v0,1.5))+xi/(rate+rate_st));
  }
  virtual inline Real getSteadyTangent(Real rate, Real state, Real a, Real b, Real D,
				       Real f0, Real rate_st, Real state_st){

    return (1+b*log(1+D/(sqrt(rate*rate+v0*v0)*state_st)))*(theta*v0*v0/(pow((rate*rate+v0*v0),1.5))+xi/(rate+rate_st))
      - b*D*rate/(state_st*(pow((rate*rate+v0*v0),1.5))*(1+D/(sqrt(rate*rate+v0*v0)*state_st)))*(theta/(sqrt(1+(v0/rate)*(v0/rate)))+xi*log(1+rate/rate_st));
  }
  virtual inline Real getStableState(Real interface_traction, Real sigma_0, Real rate,
				     Real a, Real b, Real D, Real f0, Real rate_st, Real state_st) {
    Real ln = interface_traction/(b*sigma_0*(theta/sqrt(1+v0*v0/(rate*rate))+xi*log(1+rate/rate_st)))-1/b;
    return state_st*(exp(ln)-1);
  }
  
public:
  Real v0;
  Real theta;
  Real xi;
};
