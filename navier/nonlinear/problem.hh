#ifndef  DUNE_PROBLEM_HH__
#define  DUNE_PROBLEM_HH__


namespace Dune {

/**
 * @brief describes the initial and exact solution of the advection-diffusion model
 *
 * \f[u(x,y,z,t):=\displaystyle{\sum_{i=0}^{1}} v_i(t) \cdot \mu_i(x) \cdot
 * \nu_i(y) \cdot \omega_i(z)\f]
 *
 * with
 *
 * \f{eqnarray*}{
 * v_0(t)&:=&e^{-\varepsilon t \pi^2 (2^2 + 1^2 + 1.3^2 )} \\
 * \mu_0(x)&:=&0.6\cdot \cos(2\pi (x-at)) + 0.8\cdot \sin(2\pi (x-at)) \\
 * \nu_0(y)&:=&1.2\cdot \cos(1\pi (y-at)) + 0.4\cdot \sin(1\pi (y-at)) \\
 * v_1(t)&:=&e^{-\varepsilon t \pi^2 (0.7^2 + 0.5^2 + 0.1^2 )} \\
 * \mu_1(x)&:=&0.9\cdot \cos(0.7\pi (x-at)) + 0.2\cdot \sin(0.7\pi (x-at)) \\
 * \nu_1(y)&:=&0.3\cdot \cos(0.5\pi (y-at)) + 0.1\cdot \sin(0.5\pi (y-at))
 * \f}
 *
 * This is a solution of the AdvectionDiffusionModel for \f$g_D = u|_{\partial
 * \Omega}\f$.
 *
 */
template <class GridType>
class U0 : public Dune::Function<
            Dune::FunctionSpace<double, double,
								GridType::dimensionworld, 3>, U0<GridType> > {
 public:
  enum { ConstantVelocity = true };
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,3> RangeType;

  typedef Dune::FunctionSpace<double, double,GridType::dimensionworld, 3>
    FunctionSpaceType;
  typedef Dune::Function<FunctionSpaceType, U0<GridType> > BaseType;
  typedef double DomainFieldType;
  typedef double RangeFieldType;


  /**
   * @brief define problem parameters
   */
  U0():
    space_(),
    BaseType(space_),
    velocity_(0),
    startTime_(Parameter::getValue<double>("femhowto.startTime",0.0)),
    epsilon(Parameter::getValue<double>("femhowto.epsilon",0.1))
  {
      velocity_[0]=0.8;
      velocity_[1]=0.8;
      
      max_n_of_coefs_ = 2;
      
      //x coordinate
      common_coef_x_[0] = 2.0;
      sin_coef_x_[0] = 0.8;
      cos_coef_x_[0] = 0.6;
        
      //x coordinate
      common_coef_x_[1] = 0.7;
      sin_coef_x_[1] = 0.2;
      cos_coef_x_[1] = 0.9;
            
      //y coordinate    
      common_coef_y_[0] = 1.0;
      sin_coef_y_[0] = 0.4;
      cos_coef_y_[0] = 1.2;

    
      //y coordinate    
      common_coef_y_[1] = 0.5;
      sin_coef_y_[1] = 0.1;
      cos_coef_y_[1] = 0.3;
      
      
      //z coordinate
      common_coef_z_[0] = 1.3;
      sin_coef_z_[0] = -0.4;
      cos_coef_z_[0] = 0.1;
  
      //z coordinate    
      common_coef_z_[1] = 0.1;
      sin_coef_z_[1] = 0.2;
      cos_coef_z_[1] = -0.3;
      
      myName = "AdvectDiff";
    }


  /**
   * @brief getter for the velocity
   */
  void velocity(const DomainType& x, DomainType& v) const {
    v = velocity_;
  }
  
  /**
   * @brief evaluates \f$ u_0(x) \f$
   */
  void evaluate(const DomainType& arg, RangeType& res) const {             /*@LST0@@LST1@*/
    evaluate(arg, startTime_, res);
  }
  
  /**
   * @brief old version of the exact solution
   * 
   * old version of evaluate(const DomainType& arg, double t, RangeType& res),
   * which is still needed by the DataWriter
   */
  inline void evaluate(double t,  const DomainType& arg, RangeType& res) const { 
    evaluate(arg, t, res);
  }
    
  /**
   * @brief evaluate exact solution
   */ 
  void evaluate(const DomainType& arg, double t, RangeType& res) const {  /*@LST0@@LST1@*/
    
    res = 0.;

    for(int i=0;i<max_n_of_coefs_;++i)
      {
        if(dimDomain == 1)
          res += exp(-epsilon*t*(SQR(common_coef_x_[i]*M_PI)))\
                 *((cos_coef_x_[i]*cos(common_coef_x_[i]*M_PI*
                                       (arg[0]-velocity_[0]*t))\
                    +  sin_coef_x_[i]*sin(common_coef_x_[i]*M_PI*
                                          (arg[0]-velocity_[0]*t))));
        else if(dimDomain == 2)
          res += exp(-epsilon*t*(SQR(common_coef_x_[i]*M_PI)+
                                 SQR(common_coef_y_[i]*M_PI)))\
                 *((cos_coef_x_[i]*cos(common_coef_x_[i]*M_PI*
                                       (arg[0]-velocity_[0]*t))\
                    +  sin_coef_x_[i]*sin(common_coef_x_[i]*M_PI*
                                          (arg[0]-velocity_[0]*t)))\
                   *(cos_coef_y_[i]*cos(common_coef_y_[i]*M_PI*
                                        (arg[1]-velocity_[1]*t))\
                     + sin_coef_y_[i]*sin(common_coef_y_[i]*M_PI*
                                          (arg[1]-velocity_[1]*t))));
        else if(dimDomain == 3)
          res += exp(-epsilon*t*(SQR(common_coef_x_[i]*M_PI)+
                                 SQR(common_coef_y_[i]*M_PI)+
                                 SQR(common_coef_z_[i]*M_PI)))\
                 *((cos_coef_x_[i]*cos(common_coef_x_[i]*M_PI*
                                       (arg[0]-velocity_[0]*t))\
                    +  sin_coef_x_[i]*sin(common_coef_x_[i]*M_PI*
                                          (arg[0]-velocity_[0]*t)))\
                   *(cos_coef_y_[i]*cos(common_coef_y_[i]*M_PI*
                                        (arg[1]-velocity_[1]*t))\
                     + sin_coef_y_[i]*sin(common_coef_y_[i]*M_PI*
                                          (arg[1]-velocity_[1]*t)))\
                   *(cos_coef_z_[i]*cos(common_coef_z_[i]*M_PI*
                                        (arg[2]-velocity_[2]*t))\
                     + sin_coef_z_[i]*sin(common_coef_z_[i]*M_PI*
                                          (arg[2]-velocity_[2]*t))));
      }
  }
  
  /**
   * @brief latex output for EocOutput
   */
  std::string description()
  {
    std::ostringstream ofs;
    ofs << "Problem: " << myName << "\n\n"
      << "Epsilon = " << epsilon << "\n\n"
      << "Exact solution: $u(x,y,z,t):=\\displaystyle{\\sum_{i=0}^{" 
      << max_n_of_coefs_-1
      << "}} v_i(t) \\cdot \\mu_i(x) \\cdot \\nu_i(y) \\cdot \\omega_i(z)$\n\n";

    for(int i=0;i<max_n_of_coefs_;++i)
      { 
        std::ostringstream temp;
        if(dimDomain > 1) {
          temp << common_coef_y_[i] << "^2 ";
        }
        if(dimDomain > 2) {
          temp << "+ " << common_coef_z_[i] << "^2 ";
        }
        ofs << "$v_" << i << "(t):=e^{-\\varepsilon t \\pi^2 (" 
          << common_coef_x_[i] << "^2 + " 
          << temp.str()
          << ")} $\n\n"
          << "$\\mu_" << i << "(x):=" << cos_coef_x_[i] << "\\cdot \\cos(" << common_coef_x_[i] << "\\pi (x-at)) + " << sin_coef_x_[i] 
    << "\\cdot \\sin(" << common_coef_x_[i] << "\\pi (x-at)) $";
        if(dimDomain > 1)
        {
          ofs << "\n\n"
            << "$\\nu_" << i 
            << "(y):=" << cos_coef_y_[i] << "\\cdot \\cos(" 
            << common_coef_y_[i] << "\\pi (y-at)) + " << sin_coef_y_[i] 
            << "\\cdot \\sin(" << common_coef_y_[i] << "\\pi (y-at)) $";
        }
        if(dimDomain >2)
        {
          ofs << "\n\n"
            << "$\\omega_" << i << "(z):=" << cos_coef_z_[i] << "\\cdot \\cos(" << common_coef_z_[i] << "\\pi (z-at)) + " << sin_coef_z_[i] 
            << "\\cdot \\sin(" << common_coef_z_[i] << "\\pi (z-at)) $";
        }
        ofs << "\n\n";
      }
    ofs << "\n\n";
    return ofs.str();
  }
 private:    
  FunctionSpaceType space_;
  DomainType velocity_;
  int max_n_of_coefs_;
  double common_coef_x_[2];
  double sin_coef_x_[2];
  double cos_coef_x_[2];
  double common_coef_y_[2];
  double sin_coef_y_[2];
  double cos_coef_y_[2];
  double common_coef_z_[2];
  double sin_coef_z_[2];
  double cos_coef_z_[2];
  double startTime_;
 public:
  double epsilon;
  std::string myName;
};

};
#endif  /*DUNE_PROBLEM_HH__*/
