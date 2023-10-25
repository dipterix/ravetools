// https://github.com/mrdoob/three.js/blob/r129/src/math/Quaternion.js

#ifndef RAVETOOLS_QUATERNION_H
#define RAVETOOLS_QUATERNION_H

#include <Rcpp.h>

// #include <functional>

namespace rave3d {

class Vector3;
class Matrix4;
// class Euler;

class Quaternion {

public:
  double x;
  double y;
  double z;
  double w;

  Quaternion();

  // double operator[](unsigned int index) const;

  Quaternion& set(const double& x, const double& y, const double& z, const double& w);

  Quaternion& copy(const Quaternion& quaternion);

  // Quaternion& setFromEuler(const Euler& euler, bool update = true);
  //
  // Quaternion& setFromAxisAngle(const Vector3& axis, double angle);
  //
  // Quaternion& setFromRotationMatrix(const Matrix4& m);
  //
  // Quaternion& setFromUnitVectors(const Vector3& vFrom, const Vector3& vTo);
  //
  // [[nodiscard]] double angleTo(const Quaternion& q) const;
  //
  // Quaternion& rotateTowards(const Quaternion& q, double step);
  //
  // Quaternion& identity();
  //
  // Quaternion& invert();
  //
  // Quaternion& conjugate();
  //
  // [[nodiscard]] double dot(const Quaternion& v) const;
  //
  // [[nodiscard]] double lengthSq() const;
  //
  // [[nodiscard]] double length() const;
  //
  // Quaternion& normalize();
  //
  // Quaternion& multiply(const Quaternion& q);
  //
  // Quaternion& premultiply(const Quaternion& q);
  //
  // Quaternion& multiplyQuaternions(const Quaternion& a, const Quaternion& b);
  //
  // Quaternion& slerp(const Quaternion& qb, double t);
  //
  // [[nodiscard]] Quaternion clone() const;
  //
  // [[nodiscard]] bool equals(const Quaternion& v) const;
  //
  // bool operator==(const Quaternion& other) const;
  //
  // bool operator!=(const Quaternion& other) const;
  //
  // Quaternion& _onChange(std::function<void()> callback);
  //
  // template<class ArrayLike>
  // Quaternion& fromArray(const ArrayLike& array, unsigned int offset = 0) {
  //
  //   this->x.value_ = array[offset];
  //   this->y.value_ = array[offset + 1];
  //   this->z.value_ = array[offset + 2];
  //   this->w.value_ = array[offset + 3];
  //
  //   this->onChangeCallback_();
  //
  //   return *this;
  // }
  //
  std::vector<double> toArray();
  //
  // friend std::ostream& operator<<(std::ostream& os, const Quaternion& v) {
  //   os << "Quaternion(x=" << v.x << ", y=" << v.y << ", z=" << v.z << ", w=" << v.w << ")";
  //   return os;
  // }

// private:
//   std::function<void()> onChangeCallback_ = [] {};
};

} // namespace rave3d

#endif // RAVETOOLS_QUATERNION_H
