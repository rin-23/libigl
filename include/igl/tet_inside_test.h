// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Rinat Abdrashitov <rinat@dgp.toronto.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_POINT_INSIDE_TETS_H
#define IGL_POINT_INSIDE_TETS_H

#include "igl_inline.h"

#include <Eigen/Dense>

namespace igl
{
  // Given a point cloud determine the tets that contain each point 
  // and its corresponding barycentric coordinates
  //
  // Inputs:
  //   P  #P by 3 list of points
  //   V  #V by 3 list of tet mesh vertex positions
  //   T  #T by 4 list of tet mesh indices into rows of V
  // Output:
  //   I  #P by 1 list indiex list.  P(i) == -1 if point i is outside of the tet mesh. Otherwise
  //      its equal to the index into T
  //   B     with V so that E is also defined over S  
  template <
    typename DerivedP,
    typename DerivedV,
    typename DerivedT,
    typename DerivedI,
    typename DerivedB>
  IGL_INLINE void point_inside_tets(
    const Eigen::MatrixBase<DerivedP>& P,
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedI>& I,
    Eigen::PlainObjectBase<DerivedB>& B);
}

#ifndef IGL_STATIC_LIBRARY
#  include "tet_inside_test.cpp"
#endif

#endif
