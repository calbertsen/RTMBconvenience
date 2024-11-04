#include <TMB.h>

typedef TMBad::ad_aug ad;
typedef Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > ConstMapMatrix_NotAD;
typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > MapMatrix_NotAD;
