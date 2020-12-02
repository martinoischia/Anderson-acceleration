/*
 * test_FixedPoint.cpp
 *
 *  Created on: Feb 1, 2020
 *      Author: forma
 */
#include <cmath>
#include "FixedPointIteration.hpp"
#include "Eigen/Dense"
#include "FixedPointTraits.hpp"
#include <iostream>

int main()
  {
    using namespace apsc;
    using FixedPointIterator=FixedPointIteration<VectorTraits>;
    using IterationFunction=FixedPointIterator::IterationFunction;

    IterationFunction phi;
    // A simple iterator function that we know converges to (y(lambda), 0.739085)
    // where y(lambda)=0 if  |lambda|<1
    // if |lambda|<1 convergence to 0 , slow if |lambda| near to 1
    // at least for the non accelerated version
    // if |lambda| =1 the fixed point (0,0.739085) is unstable. In general, we do not converge.
    // if |lambda| >1 we have two fixed points, one of which is with  y=0 and is unstable, the other is stable
    //                and we converge to the second one
    // Try to change lambda and see what happens
    double lambda=0.5;
    phi=[lambda](FixedPointIterator::ArgumentType const & x){return std::vector<double>{lambda*std::sin(x[0]),std::cos(x[1])}; };
    
    FixedPointIterator iterate{phi};
    FixedPointIterator::ArgumentType startingPoint{5.0,7.0};
    std::cout<<"*** WITH BASIC METHOD:\n";
    print_result(iterate.compute(startingPoint));
    // Now with acceleration and Eigen
    using   FixedPointIterator2=apsc::FixedPointIteration<EigenTraits,ASecantAccelerator>;
    using IterationFunction2=FixedPointIterator2::IterationFunction;
    IterationFunction2 phi2;
    phi2 = [lambda](FixedPointIterator2::ArgumentType const & x)
                     {
                       FixedPointIterator2::ArgumentType res(2);
                       res[0]=lambda*std::sin(x[0]);
                       res[1]=std::cos(x[1]);
                       return res;
                     };
    FixedPointIterator2 iterate2{phi2};
    FixedPointIterator2::ArgumentType startingPoint2(2);
    startingPoint2[0]=5.0;
    startingPoint2[1]=7.0;
    std::cout<<"*** WITH SECANT ACCELERATION:\n";
    print_result(iterate2.compute(startingPoint2));
  }



