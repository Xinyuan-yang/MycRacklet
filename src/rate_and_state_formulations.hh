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
  /** Constructor. Assign the frictional parameters v0
   */
  RegularizedRandSFormulation(Real v0){
    this->v0=v0;}

  /** Standard destructor
   */
  virtual ~RegularizedRandSFormulation(){};
  /** Compute the friction coefficient 
   */
  virtual inline Real operator()(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return (1+b*log(1+state/state_st))*(f0/sqrt(1+v0*v0/(rate*rate))+a*log(1+rate/rate_st));
  }
  /** Compute and return the slope of the friction law
      @param rate - sliding velocity of a point
      @param state - state variable of a point
      @param a - R&S parameter 
      @param b - R&S parameter 
      @param D - R&S parameter 
      @param f0 - R&S parameter 
      @param rate_st - R&S parameter 
      @param state_st - R&S parameter 
   */
  virtual inline Real getTangent(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return (1+b*log(1+state/state_st))*(f0*v0*v0/(pow(rate*rate+v0*v0,1.5))+a/(rate+rate_st));
  }
  /** Compute and return the slope of the steady state friction law
   */
  virtual inline Real getSteadyTangent(Real rate, Real state, Real a, Real b, Real D,
				       Real f0, Real rate_st, Real state_st){

    return (1+b*log(1+D/(sqrt(rate*rate+v0*v0)*state_st)))*(f0*v0*v0/(pow((rate*rate+v0*v0),1.5))+a/(rate+rate_st))
      - b*D*rate/(state_st*(pow((rate*rate+v0*v0),1.5))*(1+D/(sqrt(rate*rate+v0*v0)*state_st)))*(f0/(sqrt(1+(v0/rate)*(v0/rate)))+a*log(1+rate/rate_st));
  }
  /** Compute and return the state value satisfying steady state
   */
  virtual inline Real getStableState(Real interface_traction, Real sigma_0, Real rate,
				     Real a, Real b, Real D, Real f0, Real rate_st, Real state_st) {
    Real ln = interface_traction/(b*sigma_0*(f0/sqrt(1+v0*v0/(rate*rate))+a*log(1+rate/rate_st)))-1/b;
    return state_st*(exp(ln)-1);
  }
  
public:
  /** Velocity that regularize the friction value at low velocity.
   */
  Real v0;
};

/* -------------------------------------------------------------------------- */
struct RegularizedWeakeningRandSFormulation : public RandSFormulation {
  /** Constructor. Assign the frictional parameters v0
   */
  RegularizedWeakeningRandSFormulation(Real v0){
    this->v0=v0;}

  /** Standard destructor
   */
  virtual ~RegularizedWeakeningRandSFormulation(){};
  /** Compute the friction coefficient 
   */
  virtual inline Real operator()(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return (1+b*log(state/state_st))*(f0/sqrt(1+v0*v0/(rate*rate))+a*log(1+rate/rate_st));
  }
  /** Compute and return the slope of the friction law
      @param rate - sliding velocity of a point
      @param state - state variable of a point
      @param a - R&S parameter 
      @param b - R&S parameter 
      @param D - R&S parameter 
      @param f0 - R&S parameter 
      @param rate_st - R&S parameter 
      @param state_st - R&S parameter 
   */
  virtual inline Real getTangent(Real rate, Real state, Real a, Real b, Real D,
				 Real f0, Real rate_st, Real state_st){
    return (1+b*log(state/state_st))*(f0*v0*v0/(pow(rate*rate+v0*v0,1.5))+a/(rate+rate_st));
  }
  /** Compute and return the slope of the steady state friction law
   */
  virtual inline Real getSteadyTangent(Real rate, Real state, Real a, Real b, Real D,
				       Real f0, Real rate_st, Real state_st){

    return (1+b*log(D/(sqrt(rate*rate+v0*v0)*state_st)))*(f0*v0*v0/(pow((rate*rate+v0*v0),1.5))+a/(rate+rate_st))
 - b*D*rate/(state_st*(pow((rate*rate+v0*v0),1.5))*(D/(sqrt(rate*rate+v0*v0)*state_st)))*(f0/(sqrt(1+(v0/rate)*(v0/rate)))+a*log(1+rate/rate_st));
  }
  /** Compute and return the state value satisfying steady state
   */
  virtual inline Real getStableState(Real interface_traction, Real sigma_0, Real rate,
				     Real a, Real b, Real D, Real f0, Real rate_st, Real state_st) {
    Real ln = interface_traction/(b*sigma_0*(f0/sqrt(1+v0*v0/(rate*rate))+a*log(1+rate/rate_st)))-1/b;
    return state_st*exp(ln);
  }
  
public:
  /** Velocity that regularize the friction value at low velocity.
   */
  Real v0;
};
