#include "tet_inside_test.h"
#include "barycentric_coordinates.h"


template <
  typename DerivedP>
bool insideplane(const Eigen::MatrixBase<DerivedP>& p,
                 const Eigen::MatrixBase<DerivedP>& a1,
                 const Eigen::MatrixBase<DerivedP>& a2,
                 const Eigen::MatrixBase<DerivedP>& a3) 
{   
    auto c = (a2-a1).cross(a3-a1);
    return c.dot(p-a1) <= 0;
}  

template <
  typename DerivedP,
  typename DerivedV,
  typename DerivedT,
  typename DerivedI,
  typename DerivedB>
IGL_INLINE void igl::point_inside_tets(
  const Eigen::MatrixBase<DerivedP>& P,
  const Eigen::MatrixBase<DerivedV>& V,
  const Eigen::MatrixBase<DerivedT>& T,
  Eigen::PlainObjectBase<DerivedI>& I,
  Eigen::PlainObjectBase<DerivedB>& B) 
{
    typedef Eigen::Matrix<typename DerivedV::Scalar,1,3> RowVector3S;
    assert(T.cols() == 4 && "T should be a tet mesh");
    I.resize(P.rows());
    B.resize(P.rows(), T.cols());
  
    for (int i=0; i < P.rows(); ++i) 
    {
        int inside_idx = -1;
        const RowVector3S & p = P.row(i);

        for (int j=0; j < T.rows(); ++j) 
        {
            bool inside = true;
             
            const RowVector3S & a = V.row(T(j,0));
            const RowVector3S & b = V.row(T(j,1));  
            const RowVector3S & c = V.row(T(j,2));
            const RowVector3S & d = V.row(T(j,3));
            inside &= insideplane(p, a, b, d);
            inside &= insideplane(p, a, c, b);
            inside &= insideplane(p, d, c, a);
            inside &= insideplane(p, b, c, d);

            if (inside) 
            {
                inside_idx = j;
                break; 
            }
        }

        I(i) = inside_idx; // -1 if outisde of the tet mesh
        if (inside_idx != -1) 
        {
            // RowVector3S b;
            // barycentric_coordinates(P.row(i), V.row(T(I(i),0)), V.row(T(I(i),1)), V.row(T(I(i),2)), V.row(T(I(i),3)), b);
            // B.row(i) = b; 
        }
        else 
        {
            B.row(i) << 0.0,0.0,0.0,0.0; 
        }
    }

}