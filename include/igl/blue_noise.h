// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Rinat Abdrashitov <rinat@dgp.toronto.edu>// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_BLUE_NOISE_H
#define IGL_BLUE_NOISE_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl 
{     
  // BLUENOISE generate samples in N dimensional space from a blue noise distribution.
  // based on Fast Poisson Disk Sampling in Arbitrary Dimensions by Robert Bridson
  // https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
  //
  // Parts of the code are taken from the Robert Bridson's code for curl noise
  // https://www.cs.ubc.ca/~rbridson/download/curlnoise.tar.gz
  //
  //
  // Input:
  // N dimension of the sample domain
  // radius minimum distance between samples
  // xmin, xmax extent of the sample domain
  // seed random seed
  // max_sample_attempts limit of samples to choose before rejection in the algorithm, default value is 30
  // 
  // Output:
  //  S #samples by N matrix of samples
  //
  // Example: How to sample a unit cube?
  // double radius = 0.1;
  // Vector3d xmin(-1,-1,-1), xmax(1,1,1); 
  // MatrixXd S;
  // igl::blue_noise<3>(radius, xmin, xmax, 0, 30, S);

 template <unsigned int N, typename DerivedX, typename DerivedS>
 IGL_INLINE void blue_noise(typename DerivedS::Scalar radius, 
                            const Eigen::PlainObjectBase<DerivedX>& xmin,  
                            const Eigen::PlainObjectBase<DerivedX>& xmax,  
                            unsigned int seed, 
                            int max_sample_attempts, 
                            Eigen::PlainObjectBase<DerivedS>& S);
}

#ifndef IGL_STATIC_LIBRARY
#  include "blue_noise.cpp"
#endif

#endif