// complex.hpp
//
// This file develops a class to represent complex numbers.  It is
// based (at least at present) on the STL std::complex<T> class
// template, and does not seek to provide any new functionality, but
// rather just a convenient implementation.  In particular, I have
// implemented abs(), arg(), and norm() from the std namespace as
// member functions, which is what I wish the template designers had
// done.  All member functions are defined inline inside the class
// definition, and thus, with compiler optimizations, should hopefully
// not add any significant overhead over using the raw STL template.
// Also, everything needed is in this header file - i.e., there is no
// corresponding .cpp module.
// 
#ifndef COMPLEX_H_
#define COMPLEX_H_
//

#include <complex>
#include <cmath>
#include <string>
#include "typedefs.hpp"     /* typedef Real */
#include "geom_r3.hpp"      /* XYZ, Matrix */



//////
// CLASS:   ::::  Complex  ::::
//
//   Provides a full-featured representation of complex numbers
//   consisting of a real and imaginary part utilizing the same
//   machine representation as defined by type Real (Which is presumed
//   to have been previously defined).
//
//   Implementation is as a derived class from the std::complex<T>
//   template from the STL library, which might be a horrendous idea,
//   which I will attempt to determine by trying it anyway and seeing
//   what problems might arise.  And even if all works, I'm not sure
//   what the long-term implications will be for portability,
//   platform-independence, and future-proof'ness.  Worst case
//   scenario, this class might at some point in the future need to be
//   implemented from scratch.  But I'm gonna cross my fingers and
//   leave that for another day that hopefully won't come...
//

class Complex : public std::complex<Real> {
public:
  Complex() : std::complex<Real>() {}
  Complex(Real x) : std::complex<Real>(x) {}
  Complex(Real x, Real y) : std::complex<Real>(x,y) {}
  Complex(const std::complex<Real> & other) : std::complex<Real>(other) {}
                    // (This last one allows the derived class to take
                    // assignment from objects of the base class.)

  Real abs() {return std::abs(*this);}
  Real arg() {return std::arg(*this);}
  Real norm() {return std::norm(*this);}  
                    // Norm here means magnitude-squared, which I'm
                    // not totally a fan of because it's ambiguous,
                    // but I'll stick with it because it is apparently
                    // the standard.
  
};
namespace C3 { 
  class ComplexMatrix;

  class ComplexXYZ;
}


namespace C3{    //  Definition of namespace C3

class ComplexXYZ {
  protected:
      Complex mX;
      Complex mY;
      Complex mZ;
  
  public:
      // :::::::::::::::::::::::::::::::::::::
      // ::: Constructors (ComplexXYZ Class) :::
      // :::::::::::::::::::::::::::::::::::::
  
      ComplexXYZ() : mX(0), mY(0), mZ(0) {}  // Default to zero-vector
  
      ComplexXYZ(Complex x, Complex y, Complex z)
          : mX(x), mY(y), mZ(z) {}  // Construct from complex coordinates
  
      ComplexXYZ(const R3::XYZ &realVec)  // Convert real XYZ to ComplexXYZ
          : mX(realVec.x()), mY(realVec.y()), mZ(realVec.z()) {}
  
      // :::::::::::::::::::::::::::::::::::::::
      // ::: Property-Get Methods (ComplexXYZ) ::
      // :::::::::::::::::::::::::::::::::::::::
  
      Complex x() const { return mX; }
      Complex y() const { return mY; }
      Complex z() const { return mZ; }
  
      Real MagSquared() const {
          return std::norm(mX) + std::norm(mY) + std::norm(mZ);
      }
  
      Real Mag() const {
          return std::sqrt(MagSquared());
      }
  
      Complex Dot(const ComplexXYZ &other) const {
          return (mX * std::conj(other.mX) +
                  mY * std::conj(other.mY) +
                  mZ * std::conj(other.mZ));
      }
  
      ComplexXYZ Cross(const ComplexXYZ &other) const {
          return ComplexXYZ(mY * other.mZ - mZ * other.mY,
                            mZ * other.mX - mX * other.mZ,
                            mX * other.mY - mY * other.mX);
      }
  
      ComplexXYZ ScaledBy(const Complex scale) const {
          return ComplexXYZ(scale * mX, scale * mY, scale * mZ);
      }
  
      // :::::::::::::::::::::::::::::::::::
      // ::: Operator Overloads (ComplexXYZ) ::
      // :::::::::::::::::::::::::::::::::::
  
      ComplexXYZ operator+(const ComplexXYZ &other) const {
          return ComplexXYZ(mX + other.mX, mY + other.mY, mZ + other.mZ);
      }
  
      ComplexXYZ operator-(const ComplexXYZ &other) const {
          return ComplexXYZ(mX - other.mX, mY - other.mY, mZ - other.mZ);
      }
  
      ComplexXYZ operator*(const R3::XYZ &rhs) const {
          return ComplexXYZ((mX * rhs.x()) + (mY * rhs.y()) + (mZ * rhs.z()),
                            (mX * rhs.x()) + (mY * rhs.y()) + (mZ * rhs.z()),
                            (mX * rhs.x()) + (mY * rhs.y()) + (mZ * rhs.z()));
      }
  /*
  // Function that separates the real and imaginary parts into two distinct vectors
    R3::XYZ Real() const {
      return R3::XYZ(mX.real(), mY.real(), mZ.real());
    }
    R3::XYZ Imag() const {
      return R3::XYZ(mX.imag(), mY.imag(), mZ.imag());
    }
  */    
  };




  class ComplexMatrix {
    protected:
      Complex mxx, mxy, mxz;
      Complex myx, myy, myz;
      Complex mzx, mzy, mzz;
    
    public:
      // Constructors
      ComplexMatrix(Complex xx, Complex xy, Complex xz,
                    Complex yx, Complex yy, Complex yz,
                    Complex zx, Complex zy, Complex zz) :
        mxx(xx), mxy(xy), mxz(xz),
        myx(yx), myy(yy), myz(yz),
        mzx(zx), mzy(zy), mzz(zz) {}
    
      ComplexMatrix() :
        mxx(0), mxy(0), mxz(0),
        myx(0), myy(0), myz(0),
        mzx(0), mzy(0), mzz(0) {}
    
      ~ComplexMatrix() {}
    
      // Accessors
      Complex xx() const { return mxx; }
      Complex xy() const { return mxy; }
      Complex xz() const { return mxz; }
      Complex yx() const { return myx; }
      Complex yy() const { return myy; }
      Complex yz() const { return myz; }
      Complex zx() const { return mzx; }
      Complex zy() const { return mzy; }
      Complex zz() const { return mzz; }
    
      // Transpose
      ComplexMatrix T() const {
        return ComplexMatrix(mxx, myx, mzx,
                             mxy, myy, mzy,
                             mxz, myz, mzz);
      }
    
      // Conjugate transpose (Hermitian adjoint)
      ComplexMatrix H() const {
        return ComplexMatrix(std::conj(mxx), std::conj(myx), std::conj(mzx),
                             std::conj(mxy), std::conj(myy), std::conj(mzy),
                             std::conj(mxz), std::conj(myz), std::conj(mzz));
      }
    
      // Matrix multiplication
      ComplexMatrix &operator*=(const ComplexMatrix &rhs) {
        ComplexMatrix temp;
        temp.mxx = (mxx * rhs.mxx) + (mxy * rhs.myx) + (mxz * rhs.mzx);
        temp.mxy = (mxx * rhs.mxy) + (mxy * rhs.myy) + (mxz * rhs.mzy);
        temp.mxz = (mxx * rhs.mxz) + (mxy * rhs.myz) + (mxz * rhs.mzz);
        temp.myx = (myx * rhs.mxx) + (myy * rhs.myx) + (myz * rhs.mzx);
        temp.myy = (myx * rhs.mxy) + (myy * rhs.myy) + (myz * rhs.mzy);
        temp.myz = (myx * rhs.mxz) + (myy * rhs.myz) + (myz * rhs.mzz);
        temp.mzx = (mzx * rhs.mxx) + (mzy * rhs.myx) + (mzz * rhs.mzx);
        temp.mzy = (mzx * rhs.mxy) + (mzy * rhs.myy) + (mzz * rhs.mzy);
        temp.mzz = (mzx * rhs.mxz) + (mzy * rhs.myz) + (mzz * rhs.mzz);
        *this = temp;
        return *this;
      }

      // Matrix multiplication
      ComplexMatrix &operator*=(const R3::Matrix &rhs) {
        ComplexMatrix temp;
        temp.mxx = (mxx * rhs.xx()) + (mxy * rhs.yx()) + (mxz * rhs.zx());
        temp.mxy = (mxx * rhs.xy()) + (mxy * rhs.yy()) + (mxz * rhs.zy());
        temp.mxz = (mxx * rhs.xz()) + (mxy * rhs.yz()) + (mxz * rhs.zz());
        temp.myx = (myx * rhs.xx()) + (myy * rhs.yx()) + (myz * rhs.zx());
        temp.myy = (myx * rhs.xy()) + (myy * rhs.yy()) + (myz * rhs.zy());
        temp.myz = (myx * rhs.xz()) + (myy * rhs.yz()) + (myz * rhs.zz());
        temp.mzx = (mzx * rhs.xx()) + (mzy * rhs.yx()) + (mzz * rhs.zx());
        temp.mzy = (mzx * rhs.xy()) + (mzy * rhs.yy()) + (mzz * rhs.zy());
        temp.mzz = (mzx * rhs.xz()) + (mzy * rhs.yz()) + (mzz * rhs.zz());
        *this = temp;
        return *this;
      }
      ComplexMatrix operator*(const R3::Matrix &rhs) const {
        ComplexMatrix temp;
        temp.mxx = (mxx * rhs.xx()) + (mxy * rhs.yx()) + (mxz * rhs.zx());
        temp.mxy = (mxx * rhs.xy()) + (mxy * rhs.yy()) + (mxz * rhs.zy());
        temp.mxz = (mxx * rhs.xz()) + (mxy * rhs.yz()) + (mxz * rhs.zz());
        temp.myx = (myx * rhs.xx()) + (myy * rhs.yx()) + (myz * rhs.zx());
        temp.myy = (myx * rhs.xy()) + (myy * rhs.yy()) + (myz * rhs.zy());
        temp.myz = (myx * rhs.xz()) + (myy * rhs.yz()) + (myz * rhs.zz());
        temp.mzx = (mzx * rhs.xx()) + (mzy * rhs.yx()) + (mzz * rhs.zx());
        temp.mzy = (mzx * rhs.xy()) + (mzy * rhs.yy()) + (mzz * rhs.zy());
        temp.mzz = (mzx * rhs.xz()) + (mzy * rhs.yz()) + (mzz * rhs.zz());
        return temp;
    }

      
    
      ComplexMatrix operator*(const ComplexMatrix &rhs) const {
        ComplexMatrix result = *this;
        result *= rhs;
        return result;
      }
    
      // Scalar multiplication
      ComplexMatrix &operator*=(Complex rhs) {
        mxx *= rhs; mxy *= rhs; mxz *= rhs;
        myx *= rhs; myy *= rhs; myz *= rhs;
        mzx *= rhs; mzy *= rhs; mzz *= rhs;
        return *this;
      }
    
      ComplexMatrix operator*(Complex rhs) const {
        ComplexMatrix result = *this;
        result *= rhs;
        return result;
      }
    
      // Matrix addition
      ComplexMatrix &operator+=(const ComplexMatrix &rhs) {
        mxx += rhs.mxx; mxy += rhs.mxy; mxz += rhs.mxz;
        myx += rhs.myx; myy += rhs.myy; myz += rhs.myz;
        mzx += rhs.mzx; mzy += rhs.mzy; mzz += rhs.mzz;
        return *this;
      }
    
      ComplexMatrix operator+(const ComplexMatrix &rhs) const {
        ComplexMatrix result = *this;
        result += rhs;
        return result;
      }
    
      // Magnitude (Frobenius norm)
      Real Mag() const {
        return sqrt(norm(mxx) + norm(mxy) + norm(mxz) +
                    norm(myx) + norm(myy) + norm(myz) +
                    norm(mzx) + norm(mzy) + norm(mzz));
      }

    ComplexXYZ operator*(const R3::XYZ &vec) const {
        return ComplexXYZ(
            (mxx * vec.x()) + (mxy * vec.y()) + (mxz * vec.z()),
            (myx * vec.x()) + (myy * vec.y()) + (myz * vec.z()),
            (mzx * vec.x()) + (mzy * vec.y()) + (mzz * vec.z())
        );
    }
    
    ComplexXYZ operator*(const ComplexXYZ &vec) const {
        return ComplexXYZ(
            (mxx * vec.x()) + (mxy * vec.y()) + (mxz * vec.z()),
            (myx * vec.x()) + (myy * vec.y()) + (myz * vec.z()),
            (mzx * vec.x()) + (mzy * vec.y()) + (mzz * vec.z())
        );
    }
      /*
      // Output function for debugging
      void Print() const {
        std::cout << "[" << mxx << ", " << mxy << ", " << mxz << "]\n"
                  << "[" << myx << ", " << myy << ", " << myz << "]\n"
                  << "[" << mzx << ", " << mzy << ", " << mzz << "]\n";
      }
    
      */
    }; // END CLASS ComplexMatrix
    ///
// Arithmetic overloads involving Matrices:
//
inline ComplexMatrix operator*(Complex lhs, const ComplexMatrix &rhs) {
  return rhs * lhs;    // Mult by scalar on left
}
inline ComplexMatrix operator*(Real lhs, ComplexMatrix rhs) {
  return rhs *= lhs;    // Mult by scalar on left
}
inline ComplexMatrix operator*(ComplexMatrix lhs, Real rhs) {
  return lhs *= rhs;    // Mult by scalar on right
}
inline ComplexMatrix operator+(ComplexMatrix lhs, const ComplexMatrix &rhs) {
  return lhs += rhs;    // Matrix addition
}
inline ComplexMatrix operator*(const R3::Matrix &lhs, const ComplexMatrix &rhs) {
  ComplexMatrix temp;
  temp.xx() = (lhs.xx() * rhs.xx())  + (lhs.xy() * rhs.yx()) + (lhs.xz() * rhs.zx());
  temp.xy() = (lhs.xx() * rhs.xy())  + (lhs.xy() * rhs.yy()) + (lhs.xz() * rhs.zy());
  temp.xz() = (lhs.xx() * rhs.xz())  + (lhs.xy() * rhs.yz()) + (lhs.xz() * rhs.zz());
  temp.yx() = (lhs.yx() * rhs.xx())  + (lhs.yy() * rhs.yx()) + (lhs.yz() * rhs.zx());
  temp.yy() = (lhs.yx() * rhs.xy())  + (lhs.yy() * rhs.yy()) + (lhs.yz() * rhs.zy());
  temp.yz() = (lhs.yx() * rhs.xz())  + (lhs.yy() * rhs.yz()) + (lhs.yz() * rhs.zz());
  temp.zx() = (lhs.zx() * rhs.xx())  + (lhs.zy() * rhs.yx()) + (lhs.zz() * rhs.zx());
  temp.zy() = (lhs.zx() * rhs.xy())  + (lhs.zy() * rhs.yy()) + (lhs.zz() * rhs.zy());
  temp.zz() = (lhs.zx() * rhs.xz())  + (lhs.zy() * rhs.yz()) + (lhs.zz() * rhs.zz());
  return temp;
}

}// END NAMESPACE C3

///
#endif //#infdef COMPLEX_H_
//
