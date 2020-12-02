/* -------------------------------------------------------------------------- */
struct ViscoelasticFormulation {
  ViscoelasticFormulation(){};
  virtual ~ViscoelasticFormulation(){};
  virtual inline Real getStrength(Real strength_v0, Real rate, Real vel_lim){
    return strength_v0*(1+std::max(0.,rate)/vel_lim);
  }
  virtual inline Real getTangent(Real strength_v0, Real rate, Real vel_lim){
    Real slope;
    if(rate < 0){
      slope = 0.;
    }else{
      slope = strength_v0/vel_lim;
    }
    return slope;
  }
};

/* -------------------------------------------------------------------------- */
struct ViscoelasticQuadraticFormulation : public ViscoelasticFormulation{
  ViscoelasticQuadraticFormulation(){};
  virtual ~ViscoelasticQuadraticFormulation(){};
  virtual inline Real getStrength(Real strength_v0, Real rate, Real vel_lim){
    return strength_v0*(1+(std::max(0.,rate)/vel_lim)*(std::max(0.,rate)/vel_lim));
  }
  virtual inline Real getTangent(Real strength_v0, Real rate, Real vel_lim){
    Real slope;
    if(rate < 0){
      slope = 0.;
    }else{
      slope = 2*strength_v0*rate/(vel_lim*vel_lim);
    }
    return slope;
  }
};

/* -------------------------------------------------------------------------- */
struct ViscoelasticPowerLawFormulation : public ViscoelasticFormulation {
  ViscoelasticPowerLawFormulation(){};
  virtual ~ViscoelasticPowerLawFormulation(){};
  virtual inline Real getStrength(Real strength_v0, Real rate, Real vel_lim){
    Real exp = 1 / ( 1 - std::abs(rate)/vel_lim);
    return pow(strength_v0,exp);
  }
  virtual inline Real getTangent(Real strength_v0, Real rate, Real vel_lim){
    Real slope;
    if(rate < 0){
      slope = 0.;
    }else{
      Real exp = 1 / ( 1 - std::abs(rate)/vel_lim);
      slope = vel_lim * log(strength_v0) * pow(strength_v0,exp) / pow((vel_lim - std::abs(rate)),2);
    }
    return slope;
  }
};
