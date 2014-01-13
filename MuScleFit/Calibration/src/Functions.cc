#include <MuScleFit/Calibration/interface/Functions.h>

// Service to build the scale functor corresponding to the passed identifier                                                                               
scaleFunctionBase<double * > * scaleFunctionService( const int identifier ){
  switch ( identifier ) {
  case ( 50 ): return ( new scaleFunction50<double * > ); break;
  default: std::cout << "scaleFunctionService error: wrong identifier = " << identifier << std::endl; exit(1);
  }
}


// Service to build the resolution functor corresponding to the passed identifier                                                                               
resolutionFunctionBase<double * > * resolutionFunctionService( const int identifier ){
  switch ( identifier ) {
  case ( 45 ): return ( new resolutionFunction45<double * > ); break;
  case ( 46 ): return ( new resolutionFunction46<double * > ); break;
  case ( 57 ): return ( new resolutionFunction57<double * > ); break;
  default: std::cout << "resolutionFunctService error: wrong identifier = " << identifier << std::endl; exit(1);
  }
}

