/* -------------------------------------------------------------------------- */
struct ViscoelasticFormulation {
  ViscoelasticFormulation(){};
  virtual ~ViscoelasticFormulation(){};
  virtual inline Real getStrength(Real strength_v0, Real rate, Real vel_lim){
    return strength_v0*(1+std::abs(rate)/vel_lim);
  }
  virtual inline Real getTangent(Real strength_v0, Real rate, Real vel_lim){
    return strength_v0/vel_lim;
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
    Real exp = 1 / ( 1 - std::abs(rate)/vel_lim);
    return vel_lim * log(strength_v0) * pow(strength_v0,exp) / pow((vel_lim - std::abs(rate)),2);
  }
};
