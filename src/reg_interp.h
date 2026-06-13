#ifndef RAVETOOLS_REG_INTERP_H
#define RAVETOOLS_REG_INTERP_H

// Eigen-free voxel sampling helpers shared by the resampler (resample3D.cpp)
// and the registration core (reg_core.h). Volumes are column-major (x fastest).

#include <cmath>
#include <cstddef>

namespace ravereg {

// Trilinear interpolation at continuous 0-based voxel coordinates.
// Returns false (value 0) when outside the [0, n-1] box on any axis; otherwise
// interpolates and clamps the upper neighbor at the border. When `grad` is
// non-null it receives d value / d voxel (3 components, voxel units).
//
// Defined for arithmetic element types (e.g. double, int); not for compound
// types such as Rcomplex.
template <typename T>
inline bool trilinearSample(const T* data,
                            const int nx, const int ny, const int nz,
                            const double cx, const double cy, const double cz,
                            double& value, double* grad = nullptr)
{
  if (cx < 0.0 || cy < 0.0 || cz < 0.0 ||
      cx > static_cast<double>(nx - 1) ||
      cy > static_cast<double>(ny - 1) ||
      cz > static_cast<double>(nz - 1)) {
    value = 0.0;
    if (grad) { grad[0] = grad[1] = grad[2] = 0.0; }
    return false;
  }

  const int x0 = static_cast<int>(std::floor(cx));
  const int y0 = static_cast<int>(std::floor(cy));
  const int z0 = static_cast<int>(std::floor(cz));
  const int x1 = (x0 < nx - 1) ? x0 + 1 : x0;
  const int y1 = (y0 < ny - 1) ? y0 + 1 : y0;
  const int z1 = (z0 < nz - 1) ? z0 + 1 : z0;

  const double fx = cx - x0;
  const double fy = cy - y0;
  const double fz = cz - z0;

  const std::size_t sy = static_cast<std::size_t>(nx);
  const std::size_t sz = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);

  const double c000 = static_cast<double>(data[x0 + sy * y0 + sz * z0]);
  const double c100 = static_cast<double>(data[x1 + sy * y0 + sz * z0]);
  const double c010 = static_cast<double>(data[x0 + sy * y1 + sz * z0]);
  const double c110 = static_cast<double>(data[x1 + sy * y1 + sz * z0]);
  const double c001 = static_cast<double>(data[x0 + sy * y0 + sz * z1]);
  const double c101 = static_cast<double>(data[x1 + sy * y0 + sz * z1]);
  const double c011 = static_cast<double>(data[x0 + sy * y1 + sz * z1]);
  const double c111 = static_cast<double>(data[x1 + sy * y1 + sz * z1]);

  const double c00 = c000 + fx * (c100 - c000);
  const double c10 = c010 + fx * (c110 - c010);
  const double c01 = c001 + fx * (c101 - c001);
  const double c11 = c011 + fx * (c111 - c011);
  const double c0 = c00 + fy * (c10 - c00);
  const double c1 = c01 + fy * (c11 - c01);
  value = c0 + fz * (c1 - c0);

  if (grad) {
    const double d00 = c100 - c000;
    const double d10 = c110 - c010;
    const double d01 = c101 - c001;
    const double d11 = c111 - c011;
    const double dx0 = d00 + fy * (d10 - d00);
    const double dx1 = d01 + fy * (d11 - d01);
    grad[0] = dx0 + fz * (dx1 - dx0);
    const double ey0 = c10 - c00;
    const double ey1 = c11 - c01;
    grad[1] = ey0 + fz * (ey1 - ey0);
    grad[2] = c1 - c0;
  }
  return true;
}

// Keys cubic (Catmull-Rom, alpha=-0.5) interpolation at continuous 0-based voxel
// coordinates. Out-of-bounds (outside [0,N-1] on any axis) returns false.
// Neighbor indices outside bounds are clamped (zero-flux Neumann / edge extension).
// When `grad` is non-null it receives d value / d voxel (3 components, voxel units).
template <typename T>
inline bool bsplineSample(const T* data,
                          const int nx, const int ny, const int nz,
                          const double cx, const double cy, const double cz,
                          double& value, double* grad = nullptr)
{
  if (cx < 0.0 || cy < 0.0 || cz < 0.0 ||
      cx > static_cast<double>(nx - 1) ||
      cy > static_cast<double>(ny - 1) ||
      cz > static_cast<double>(nz - 1)) {
    value = 0.0;
    if (grad) grad[0] = grad[1] = grad[2] = 0.0;
    return false;
  }

  const int x0 = static_cast<int>(std::floor(cx));
  const int y0 = static_cast<int>(std::floor(cy));
  const int z0 = static_cast<int>(std::floor(cz));

  // 4-tap Catmull-Rom weights for each axis.
  // For fractional offset f the four distances are (1+f, f, 1-f, 2-f).
  // Branch [0,1]: 1.5d³-2.5d²+1; branch [1,2]: -0.5d³+2.5d²-4d+2.
  double wx[4], wy[4], wz[4];
  double dwx[4], dwy[4], dwz[4];

  {
    const double f=cx-x0, d0=1.0+f, d1=f, d2=1.0-f, d3=2.0-f;
    wx[0]=d0*d0*(-0.5*d0+2.5)-4.0*d0+2.0;
    wx[1]=d1*d1*(1.5*d1-2.5)+1.0;
    wx[2]=d2*d2*(1.5*d2-2.5)+1.0;
    wx[3]=d3*d3*(-0.5*d3+2.5)-4.0*d3+2.0;
    if (grad) {
      dwx[0]=d0*(-1.5*d0+5.0)-4.0; dwx[1]=d1*(4.5*d1-5.0);
      dwx[2]=d2*(-4.5*d2+5.0);     dwx[3]=d3*(1.5*d3-5.0)+4.0;
    }
  }
  {
    const double f=cy-y0, d0=1.0+f, d1=f, d2=1.0-f, d3=2.0-f;
    wy[0]=d0*d0*(-0.5*d0+2.5)-4.0*d0+2.0;
    wy[1]=d1*d1*(1.5*d1-2.5)+1.0;
    wy[2]=d2*d2*(1.5*d2-2.5)+1.0;
    wy[3]=d3*d3*(-0.5*d3+2.5)-4.0*d3+2.0;
    if (grad) {
      dwy[0]=d0*(-1.5*d0+5.0)-4.0; dwy[1]=d1*(4.5*d1-5.0);
      dwy[2]=d2*(-4.5*d2+5.0);     dwy[3]=d3*(1.5*d3-5.0)+4.0;
    }
  }
  {
    const double f=cz-z0, d0=1.0+f, d1=f, d2=1.0-f, d3=2.0-f;
    wz[0]=d0*d0*(-0.5*d0+2.5)-4.0*d0+2.0;
    wz[1]=d1*d1*(1.5*d1-2.5)+1.0;
    wz[2]=d2*d2*(1.5*d2-2.5)+1.0;
    wz[3]=d3*d3*(-0.5*d3+2.5)-4.0*d3+2.0;
    if (grad) {
      dwz[0]=d0*(-1.5*d0+5.0)-4.0; dwz[1]=d1*(4.5*d1-5.0);
      dwz[2]=d2*(-4.5*d2+5.0);     dwz[3]=d3*(1.5*d3-5.0)+4.0;
    }
  }

  int ix[4], iy[4], iz[4];
  ix[0]=std::max(0,x0-1); ix[1]=x0;
  ix[2]=std::min(nx-1,x0+1); ix[3]=std::min(nx-1,x0+2);
  iy[0]=std::max(0,y0-1); iy[1]=y0;
  iy[2]=std::min(ny-1,y0+1); iy[3]=std::min(ny-1,y0+2);
  iz[0]=std::max(0,z0-1); iz[1]=z0;
  iz[2]=std::min(nz-1,z0+1); iz[3]=std::min(nz-1,z0+2);

  const std::size_t sy = static_cast<std::size_t>(nx);
  const std::size_t sz = static_cast<std::size_t>(nx)*static_cast<std::size_t>(ny);

  double v=0.0, gx=0.0, gy=0.0, gz=0.0;
  const bool need_grad = (grad != nullptr);
  for (int kk=0; kk<4; ++kk) {
    for (int jj=0; jj<4; ++jj) {
      const std::size_t syjk = sy*iy[jj] + sz*iz[kk];
      for (int ii=0; ii<4; ++ii) {
        const double d = static_cast<double>(data[ix[ii]+syjk]);
        v += wx[ii]*wy[jj]*wz[kk]*d;
        if (need_grad) {
          gx += dwx[ii]*wy[jj]*wz[kk]*d;
          gy += wx[ii]*dwy[jj]*wz[kk]*d;
          gz += wx[ii]*wy[jj]*dwz[kk]*d;
        }
      }
    }
  }
  value = v;
  if (need_grad) { grad[0]=gx; grad[1]=gy; grad[2]=gz; }
  return true;
}

} // namespace ravereg

#endif // RAVETOOLS_REG_INTERP_H
