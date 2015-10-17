#ifndef INC_Matrix33_H
#define INC_Matrix33_H

#include "MatrixNM.h"
#include "Vector3.h"
#include "JPMathWrite.h"
#include "eig3.h"

class Matrix22;
class Matrix44;

class Matrix33 {
public:
  inline Matrix33();
  inline Matrix33(double x);
  inline Matrix33(const double &diag00, const double &diag11, const double &diag22);
  inline Matrix33(const Vector3 &diag);
  inline Matrix33(const Matrix33 &that);
  Matrix33(const Matrix22 &that);
  Matrix33(const Matrix44 &that);
  Matrix33(const MatrixNM &that);
  inline Matrix33(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3);

  inline Matrix33(double M00, double M01, double M02,
				  double M10, double M11, double M12,
				  double M20, double M21, double M22);

  inline Matrix33(double M[3][3]);

  inline Matrix33 &set      (const double d);
  inline Matrix33 &operator=(const double d);

  inline Matrix33 &set      (const Matrix33 &that);  
  inline Matrix33 &operator=(const Matrix33 &that); 

	inline Matrix33 &set			 (double M[3][3]);
  inline Matrix33 &operator=(double M[3][3]); 

	inline Matrix33 &set     (const Vector3 &v1, const Vector3 &v2, const Vector3 &v3);
  inline Matrix33 &set			(double M00, double M01, double M02,
					      							 double M10, double M11, double M12,
															 double M20, double M21, double M22);
  inline void setColumn(int j, Vector3 v);

  inline int operator!=(const Matrix33 &that)const; 
  inline int operator==(const Matrix33 &that)const; 

  inline int operator==(double d) const;
  inline int operator!=(double d) const;
  
  inline Matrix33 &operator+=(double d);
  inline Matrix33 &operator-=(double d);
  inline Matrix33 &operator*=(double d);
  inline Matrix33 &operator/=(double d);

  // component-wise operations.
  inline Matrix33 &operator+=(const Matrix33 &that);
  inline Matrix33 &operator-=(const Matrix33 &that);
  //inline Matrix33 &operator*=(const Matrix33 &that);

  // Left : this = that x this  
  // Right: this = this x that
  Matrix33 &leftMultiply (const Matrix33 &that);
  Matrix33 &rightMultiply(const Matrix33 &that);

  inline Matrix33 &setIdentity     ();

  inline void writeMatrix();

  inline void eigVecVals(double *e1, double *e2, double *e3, Vector3 *ev1, Vector3 *ev2, Vector3 *ev3);
  inline void eigVecVals(double eval[3], Vector3 evec[3]);

public:
  union{
	  double elements[3][3];
	  struct{
		  double m00, m01, m02, m10, m11, m12, m20, m21, m22;
	  };
  };

};

inline Matrix33 operator-(const Matrix33 &v);

inline Matrix33 operator+(const Matrix33 &a, double b);
inline Matrix33 operator-(const Matrix33 &a, double b);
inline Matrix33 operator*(const Matrix33 &a, double b);

inline Matrix33 operator+(const Matrix33 &a, const Matrix33 &b);
inline Matrix33 operator-(const Matrix33 &a, const Matrix33 &b);
inline Matrix33 operator*(const Matrix33 &a, const Matrix33 &b); 

inline Matrix33 multiply(const Matrix33 &a, const Matrix33 &b); 
inline Matrix33 componenentMult(const Matrix33 &a, const Matrix33 &b); 
inline Matrix33 conjugate(const Matrix33 &a, const Matrix33 &b); 
inline Matrix33 othoconjugate(const Matrix33 &a, const Matrix33 &b); 
inline Vector3   operator*(const Matrix33 &a, const Vector3   &b);
inline Vector3   operator*(const Vector3   &a, const Matrix33 &b);

inline double determinant(const Matrix33 &a);

inline Matrix33 transpose(const Matrix33 &a);
inline Matrix33 symmetricAdjoint(Matrix33 &M);
inline Matrix33   inverse(const Matrix33 &a);
inline Matrix33 outer(const Vector3 &v1, const Vector3 &v2);
inline Matrix33 rotMat3(Vector3 v, double angle);

inline void normalize_F(Matrix33 &a);

inline double l_frobenius(Matrix33 &m);
inline double l_infinity(Matrix33 &m);

inline Matrix33::Matrix33() {
  elements[0][0] = 1;
  elements[0][1] = 0;
  elements[0][2] = 0;
  elements[1][0] = 0;
  elements[1][1] = 1;
  elements[1][2] = 0;
  elements[2][0] = 0;
  elements[2][1] = 0;
  elements[2][2] = 1;
}

inline Matrix33::Matrix33(double x) {
  elements[0][0] = x;
  elements[0][1] = x;
  elements[0][2] = x;
  elements[1][0] = x;
  elements[1][1] = x;
  elements[1][2] = x;
  elements[2][0] = x;
  elements[2][1] = x;
  elements[2][2] = x;
}

inline Matrix33::Matrix33(const double &diag00, const double &diag11, const double &diag22){
	elements[0][0] = diag00;
	elements[0][1] = 0;
	elements[0][2] = 0;
	elements[1][0] = 0;
	elements[1][1] = diag11;
	elements[1][2] = 0;
	elements[2][0] = 0;
	elements[2][1] = 0;
	elements[2][2] = diag22;
}

inline Matrix33::Matrix33(const Vector3 &diag){
	elements[0][0] = diag.a;
	elements[0][1] = 0;
	elements[0][2] = 0;
	elements[1][0] = 0;
	elements[1][1] = diag.b;
	elements[1][2] = 0;
	elements[2][0] = 0;
	elements[2][1] = 0;
	elements[2][2] = diag.c;
}

inline Matrix33::Matrix33(double M00, double M01, double M02,
																double M10, double M11, double M12,
																double M20, double M21, double M22) {
  elements[0][0] = M00;
  elements[0][1] = M01;
  elements[0][2] = M02;
  elements[1][0] = M10;
  elements[1][1] = M11;
  elements[1][2] = M12;
  elements[2][0] = M20;
  elements[2][1] = M21;
  elements[2][2] = M22;
};

inline Matrix33::Matrix33(const Matrix33 &that) {
  elements[0][0] = that.elements[0][0];
  elements[0][1] = that.elements[0][1];
  elements[0][2] = that.elements[0][2];
  elements[1][0] = that.elements[1][0];
  elements[1][1] = that.elements[1][1];
  elements[1][2] = that.elements[1][2];
  elements[2][0] = that.elements[2][0];
  elements[2][1] = that.elements[2][1];
  elements[2][2] = that.elements[2][2];
};

inline Matrix33::Matrix33(const MatrixNM &that){
	if(that.n == 3 && that.m == 3){
		elements[0][0] = that.elements[0][0];
		elements[0][1] = that.elements[0][1];
		elements[0][2] = that.elements[0][2];
		elements[1][0] = that.elements[1][0];
		elements[1][1] = that.elements[1][1];
		elements[1][2] = that.elements[1][2];
		elements[2][0] = that.elements[2][0];
		elements[2][1] = that.elements[2][1];
		elements[2][2] = that.elements[2][2];
	}
};

inline Matrix33::Matrix33(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3) {
	/*elements[0][0] = v1.a;
	elements[0][1] = v1.b;
	elements[0][2] = v1.c;
	elements[1][0] = v2.a;
	elements[1][1] = v2.b;
	elements[1][2] = v2.c;
	elements[2][0] = v3.a;
	elements[2][1] = v3.b;
	elements[2][2] = v3.c;*/

	elements[0][0] = v1.a;
	elements[1][0] = v1.b;
	elements[2][0] = v1.c;
	elements[0][1] = v2.a;
	elements[1][1] = v2.b;
	elements[2][1] = v2.c;
	elements[0][2] = v3.a;
	elements[1][2] = v3.b;
	elements[2][2] = v3.c;
}

inline Matrix33 &Matrix33::set(const double d) {
  return (*this)=d;
}

inline Matrix33 &Matrix33::operator=(const double d) {
  elements[0][0] = d;
  elements[0][1] = d;
  elements[0][2] = d;

  elements[1][0] = d;
  elements[1][1] = d;
  elements[1][2] = d;

  elements[2][0] = d;
  elements[2][1] = d;
  elements[2][2] = d;

  return (*this);
};

inline void Matrix33::setColumn(int j, Vector3 v){
	elements[0][j] = v.elements[0];
	elements[1][j] = v.elements[1];
	elements[2][j] = v.elements[2];
}

inline Matrix33 &Matrix33::set(const Matrix33 &that) {
  return (*this)=that;
}

inline Matrix33 &Matrix33::operator=(const Matrix33 &that) {
  elements[0][0] = that.elements[0][0];
  elements[0][1] = that.elements[0][1];
  elements[0][2] = that.elements[0][2];
  elements[1][0] = that.elements[1][0];
  elements[1][1] = that.elements[1][1];
  elements[1][2] = that.elements[1][2];
  elements[2][0] = that.elements[2][0];
  elements[2][1] = that.elements[2][1];
  elements[2][2] = that.elements[2][2];
  return (*this);
};

inline Matrix33 &Matrix33::set(double M[3][3]) {
  return (*this)=M;
}

inline Matrix33 &Matrix33::operator=(double M[3][3]) {
  elements[0][0] = M[0][0];
  elements[0][1] = M[0][1];
  elements[0][2] = M[0][2];

  elements[1][0] = M[1][0];
  elements[1][1] = M[1][1];
  elements[1][2] = M[1][2];

  elements[2][0] = M[2][0];
  elements[2][1] = M[2][1];
  elements[2][2] = M[2][2];
return (*this);
};

inline Matrix33 &Matrix33::set(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3) {
	/*elements[0][0] = v1.a;
	elements[0][1] = v1.b;
	elements[0][2] = v1.c;
	elements[1][0] = v2.a;
	elements[1][1] = v2.b;
	elements[1][2] = v2.c;
	elements[2][0] = v3.a;
	elements[2][1] = v3.b;
	elements[2][2] = v3.c;*/

	elements[0][0] = v1.a;
	elements[1][0] = v1.b;
	elements[2][0] = v1.c;
	elements[0][1] = v2.a;
	elements[1][1] = v2.b;
	elements[2][1] = v2.c;
	elements[0][2] = v3.a;
	elements[1][2] = v3.b;
	elements[2][2] = v3.c;
	return (*this);
}

inline Matrix33 &Matrix33::set			(double M00, double M01, double M02,
				      							 double M10, double M11, double M12,
														 double M20, double M21, double M22)
{
	elements[0][0] = M00;
	elements[0][1] = M01;
	elements[0][2] = M02;
	elements[1][0] = M10;
	elements[1][1] = M11;
	elements[1][2] = M12;
	elements[2][0] = M20;
	elements[2][1] = M21;
	elements[2][2] = M22;
	return (*this);
}

inline int Matrix33::operator==(double d) const {
  return  ( (elements[0][0] == d) && (elements[0][1] == d) && (elements[0][2] == d) &&
						(elements[1][0] == d) && (elements[1][1] == d) && (elements[1][2] == d) && 
						(elements[2][0] == d) && (elements[2][1] == d) && (elements[2][2] == d));
}

inline int Matrix33::operator!=(double d) const {
  return  ( (elements[0][0] != d) || (elements[0][1] != d) || (elements[0][2] != d) ||
						(elements[1][0] != d) || (elements[1][1] != d) || (elements[1][2] != d) ||
						(elements[2][0] != d) || (elements[2][1] != d) || (elements[2][2] != d));
}
  
inline int Matrix33::operator==(const Matrix33 &that)const {
  return ( (elements[0][0] == that.elements[0][0]) && (elements[0][1] == that.elements[0][1]) && (elements[0][2] == that.elements[0][2]) &&
					 (elements[1][0] == that.elements[1][0]) && (elements[1][1] == that.elements[1][1]) && (elements[1][2] == that.elements[1][2]) &&
					 (elements[2][0] == that.elements[2][0]) && (elements[2][1] == that.elements[2][1]) && (elements[2][2] == that.elements[2][2]));
}

inline int Matrix33::operator!=(const Matrix33 &that)const {
  return ( (elements[0][0] != that.elements[0][0]) || (elements[0][1] != that.elements[0][1]) || (elements[0][2] != that.elements[0][2]) ||
					 (elements[1][0] != that.elements[1][0]) || (elements[1][1] != that.elements[1][1]) || (elements[1][2] != that.elements[1][2]) ||
					 (elements[2][0] != that.elements[2][0]) || (elements[2][1] != that.elements[2][1]) || (elements[2][2] != that.elements[2][2]));
}

inline Matrix33 &Matrix33::operator+=(double d) {
  elements[0][0] += d; elements[0][1] += d; elements[0][2] += d; 
  elements[1][0] += d; elements[1][1] += d; elements[1][2] += d; 
  elements[2][0] += d; elements[2][1] += d; elements[2][2] += d; 
  return (*this);
}

inline Matrix33 &Matrix33::operator-=(double d) {
  elements[0][0] -= d; elements[0][1] -= d; elements[0][2] -= d; 
  elements[1][0] -= d; elements[1][1] -= d; elements[1][2] -= d;
  elements[2][0] -= d; elements[2][1] -= d; elements[2][2] -= d;
  return (*this);
}

inline Matrix33 &Matrix33::operator*=(double d) {
  elements[0][0] *= d; elements[0][1] *= d; elements[0][2] *= d; 
  elements[1][0] *= d; elements[1][1] *= d; elements[1][2] *= d; 
  elements[2][0] *= d; elements[2][1] *= d; elements[2][2] *= d; 
  return (*this);
}

inline Matrix33 &Matrix33::operator/=(double d) {
  elements[0][0] /= d; elements[0][1] /= d; elements[0][2] /= d; 
  elements[1][0] /= d; elements[1][1] /= d; elements[1][2] /= d; 
  elements[2][0] /= d; elements[2][1] /= d; elements[2][2] /= d; 
  return (*this);
}

inline Matrix33 &Matrix33::operator+=(const Matrix33 &that) {
  elements[0][0] += that.elements[0][0]; elements[0][1] += that.elements[0][1]; elements[0][2] += that.elements[0][2]; 
  elements[1][0] += that.elements[1][0]; elements[1][1] += that.elements[1][1]; elements[1][2] += that.elements[1][2]; 
  elements[2][0] += that.elements[2][0]; elements[2][1] += that.elements[2][1]; elements[2][2] += that.elements[2][2]; 
  return (*this);
}
  
inline Matrix33 &Matrix33::operator-=(const Matrix33 &that) {
  elements[0][0] -= that.elements[0][0]; elements[0][1] -= that.elements[0][1]; elements[0][2] -= that.elements[0][2]; 
  elements[1][0] -= that.elements[1][0]; elements[1][1] -= that.elements[1][1]; elements[1][2] -= that.elements[1][2]; 
  elements[2][0] -= that.elements[2][0]; elements[2][1] -= that.elements[2][1]; elements[2][2] -= that.elements[2][2]; 
  return (*this);
}

/*inline Matrix33 &Matrix33::operator*=(const Matrix33 &that) {
  elements[0][0] *= that.elements[0][0]; elements[0][1] *= that.elements[0][1]; elements[0][2] *= that.elements[0][2]; 
  elements[1][0] *= that.elements[1][0]; elements[1][1] *= that.elements[1][1]; elements[1][2] *= that.elements[1][2]; 
  elements[2][0] *= that.elements[2][0]; elements[2][1] *= that.elements[2][1]; elements[2][2] *= that.elements[2][2]; 
  return (*this);
}*/

inline Matrix33 &Matrix33::leftMultiply (const Matrix33 &that){
	Matrix33 tmp(elements[0][0], elements[0][1], elements[0][2], 
									elements[1][0], elements[1][1], elements[1][2],
									elements[2][0], elements[2][1], elements[2][2]);
	
	elements[0][0] = that.elements[0][0] * tmp.elements[0][0] + that.elements[0][1] * tmp.elements[1][0] + that.elements[0][2] * tmp.elements[2][0];
	elements[0][1] = that.elements[0][0] * tmp.elements[0][1] + that.elements[0][1] * tmp.elements[1][1] + that.elements[0][2] * tmp.elements[2][1];
	elements[0][2] = that.elements[0][0] * tmp.elements[0][2] + that.elements[0][1] * tmp.elements[1][2] + that.elements[0][2] * tmp.elements[2][2];

	elements[1][0] = that.elements[1][0] * tmp.elements[0][0] + that.elements[1][1] * tmp.elements[1][0] + that.elements[1][2] * tmp.elements[2][0];
	elements[1][1] = that.elements[1][0] * tmp.elements[0][1] + that.elements[1][1] * tmp.elements[1][1] + that.elements[1][2] * tmp.elements[2][1];
	elements[1][2] = that.elements[1][0] * tmp.elements[0][2] + that.elements[1][1] * tmp.elements[1][2] + that.elements[1][2] * tmp.elements[2][2];

	elements[2][0] = that.elements[2][0] * tmp.elements[0][0] + that.elements[2][1] * tmp.elements[1][0] + that.elements[2][2] * tmp.elements[2][0];
	elements[2][1] = that.elements[2][0] * tmp.elements[0][1] + that.elements[2][1] * tmp.elements[1][1] + that.elements[2][2] * tmp.elements[2][1];
	elements[2][2] = that.elements[2][0] * tmp.elements[0][2] + that.elements[2][1] * tmp.elements[1][2] + that.elements[2][2] * tmp.elements[2][2];
	return (*this);
};

inline Matrix33 &Matrix33::rightMultiply(const Matrix33 &that){
	Matrix33 tmp(elements[0][0], elements[0][1], elements[0][2], 
									elements[1][0], elements[1][1], elements[1][2],
									elements[2][0], elements[2][1], elements[2][2]);

	elements[0][0] = tmp.elements[0][0] * that.elements[0][0] + tmp.elements[0][1] * that.elements[1][0] + tmp.elements[0][2] * that.elements[2][0];
	elements[0][1] = tmp.elements[0][0] * that.elements[0][1] + tmp.elements[0][1] * that.elements[1][1] + tmp.elements[0][2] * that.elements[2][1];
	elements[0][2] = tmp.elements[0][0] * that.elements[0][2] + tmp.elements[0][1] * that.elements[1][2] + tmp.elements[0][2] * that.elements[2][2];

	elements[1][0] = tmp.elements[1][0] * that.elements[0][0] + tmp.elements[1][1] * that.elements[1][0] + tmp.elements[1][2] * that.elements[2][0];
	elements[1][1] = tmp.elements[1][0] * that.elements[0][1] + tmp.elements[1][1] * that.elements[1][1] + tmp.elements[1][2] * that.elements[2][1];
	elements[1][2] = tmp.elements[1][0] * that.elements[0][2] + tmp.elements[1][1] * that.elements[1][2] + tmp.elements[1][2] * that.elements[2][2];

	elements[2][0] = tmp.elements[2][0] * that.elements[0][0] + tmp.elements[2][1] * that.elements[1][0] + tmp.elements[2][2] * that.elements[2][0];
	elements[2][1] = tmp.elements[2][0] * that.elements[0][1] + tmp.elements[2][1] * that.elements[1][1] + tmp.elements[2][2] * that.elements[2][1];
	elements[2][2] = tmp.elements[2][0] * that.elements[0][2] + tmp.elements[2][1] * that.elements[1][2] + tmp.elements[2][2] * that.elements[2][2];
	return (*this);
};

inline Matrix33 &Matrix33::setIdentity() {
  elements[0][0] = 1; elements[0][1] = 0; elements[0][2] = 0; 
  elements[1][0] = 0; elements[1][1] = 1; elements[1][2] = 0; 
  elements[2][0] = 0; elements[2][1] = 0; elements[2][2] = 1; 
  return (*this);
};

void Matrix33::writeMatrix(){
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			JPMath_write("%10f", elements[i][j]);
		}
		JPMath_write("\n");
	}
}


inline void Matrix33::eigVecVals(double *e1, double *e2, double *e3, Vector3 *ev1, Vector3 *ev2, Vector3 *ev3){
	double eVec[3][3];
	double eVal[3];

	eigen_decomposition(elements, eVec, eVal);

	*e1 = eVal[2];
	*e2 = eVal[1];
	*e3 = eVal[0];

	ev1->set(eVec[0][2], eVec[1][2], eVec[2][2]);
	ev2->set(eVec[0][1], eVec[1][1], eVec[2][1]);
	ev3->set(eVec[0][0], eVec[1][0], eVec[2][0]);
}

inline void Matrix33::eigVecVals(double eval[3], Vector3 evec[3]){
	eigVecVals(&(eval[0]), &(eval[1]), &(eval[2]), &(evec[0]), &(evec[1]), &(evec[2]));
}


inline Matrix33 operator-(const Matrix33 &a){
	 return Matrix33(-a.elements[0][0],-a.elements[0][1],-a.elements[0][2],
					 -a.elements[1][0],-a.elements[1][1],-a.elements[1][2],
					 -a.elements[2][0],-a.elements[2][1],-a.elements[2][2]);
}

inline void normalize_F(Matrix33 &a){
	double mag = l_frobenius(a);
	if(mag > 100*EPS) a /= mag;
}

inline Matrix33 operator+(const Matrix33 &a,double b) {
  return (Matrix33(a)+=b);
}

inline Matrix33 operator-(const Matrix33 &a,double b) {
  return (Matrix33(a)-=b);
}

inline Matrix33 operator*(const Matrix33 &a,double b) {
  return  Matrix33(a)*=b;
}
 
inline Matrix33 operator+(double a, const Matrix33 &b) {
return b+a;
}

inline Matrix33 operator-(double a, const Matrix33 &b) {
  return Matrix33(a-b.elements[0][0],a-b.elements[0][1],a-b.elements[0][2],
										 a-b.elements[1][0],a-b.elements[1][1],a-b.elements[1][2],
										 a-b.elements[2][0],a-b.elements[2][1],a-b.elements[2][2]);
}

inline Matrix33 operator*(double a, const Matrix33 &b) {
  return b*a;
}
 
inline Matrix33 operator+(const Matrix33 &a,const Matrix33 &b) {
  return (Matrix33(a)+=b);
}
 
inline Matrix33 operator-(const Matrix33 &a,const Matrix33 &b) {
  return (Matrix33(a)-=b);
}

inline Matrix33 operator*(const Matrix33 &a,const Matrix33 &b) {
  return  multiply(a,b);
}

inline Matrix33 multiply(const Matrix33 &a,const Matrix33 &b) {
  Matrix33 tmp(a);
  tmp.rightMultiply(b);
  return tmp;
}

inline Matrix33 conjugate(const Matrix33 a, const Matrix33 &b) {
  Matrix33 tmp(a);
	Matrix33 c = inverse(b);
  tmp.rightMultiply(b);
	tmp.leftMultiply(c);
  return tmp;
}

inline Matrix33 othoconjugate(const Matrix33 &a, const Matrix33 &b) {
  Matrix33 tmp(a);
	Matrix33 c = transpose(b);
  tmp.rightMultiply(b);
	tmp.leftMultiply(c);
  return tmp;
}

inline Vector3 operator*(const Matrix33 &a,const Vector3 &b) {
  return Vector3(b.a*a.elements[0][0] + b.b*a.elements[0][1] + b.c*a.elements[0][2], 
									 b.a*a.elements[1][0] + b.b*a.elements[1][1] + b.c*a.elements[1][2],
									 b.a*a.elements[2][0] + b.b*a.elements[2][1] + b.c*a.elements[2][2]);
}

inline Vector3 operator*(const Vector3 &a,const Matrix33 &b) {
  return Vector3(a.a*b.elements[0][0] + a.b*b.elements[1][0] + a.c*b.elements[2][0],
									 a.a*b.elements[0][1] + a.b*b.elements[1][1] + a.c*b.elements[2][1],
									 a.a*b.elements[0][2] + a.b*b.elements[1][2] + a.c*b.elements[2][2]);
}

inline double determinant(const Matrix33 &a) {
  return ( a.elements[0][0] * a.elements[1][1] * a.elements[2][2] - a.elements[2][0] * a.elements[1][1] * a.elements[0][2]
		     + a.elements[1][0] * a.elements[2][1] * a.elements[0][2] - a.elements[0][0] * a.elements[2][1] * a.elements[1][2]
				 + a.elements[2][0] * a.elements[0][1] * a.elements[1][2] - a.elements[1][0] * a.elements[0][1] * a.elements[2][2]);
}

inline Matrix33 transpose(const Matrix33 &a) {
  Matrix33 tmp(a);

	tmp.elements[0][1] = a.elements[1][0];
	tmp.elements[1][0] = a.elements[0][1];

	tmp.elements[0][2] = a.elements[2][0];
	tmp.elements[2][0] = a.elements[0][2];

	tmp.elements[2][1] = a.elements[1][2];
	tmp.elements[1][2] = a.elements[2][1];
  return tmp;
}

inline Matrix33 inverse(const Matrix33 &a) {
	Matrix33 tmp;
	double dmt;
	
	if ((dmt=determinant(a))!= 0.0) {
		tmp.elements[0][0] = (a.elements[1][1] * a.elements[2][2] - a.elements[2][1] * a.elements[1][2])/dmt;
		tmp.elements[0][1] = (a.elements[2][1] * a.elements[0][2] - a.elements[0][1] * a.elements[2][2])/dmt;
		tmp.elements[0][2] = (a.elements[0][1] * a.elements[1][2] - a.elements[1][1] * a.elements[0][2])/dmt;

		tmp.elements[1][0] = (a.elements[1][2] * a.elements[2][0] - a.elements[2][2] * a.elements[1][0])/dmt;
		tmp.elements[1][1] = (a.elements[2][2] * a.elements[0][0] - a.elements[0][2] * a.elements[2][0])/dmt;
		tmp.elements[1][2] = (a.elements[0][2] * a.elements[1][0] - a.elements[1][2] * a.elements[0][0])/dmt;

		tmp.elements[2][0] = (a.elements[1][0] * a.elements[2][1] - a.elements[2][0] * a.elements[1][1])/dmt;
		tmp.elements[2][1] = (a.elements[2][0] * a.elements[0][1] - a.elements[0][0] * a.elements[2][1])/dmt;
		tmp.elements[2][2] = (a.elements[0][0] * a.elements[1][1] - a.elements[1][0] * a.elements[0][1])/dmt;
	}
	return tmp;
}

Matrix33 outer(const Vector3 &v1, const Vector3 &v2){
	Matrix33 op;

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			op.elements[i][j] = v1.elements[i]*v2.elements[j];
		}
	}
	return op;
}

inline Matrix33 pow(const Matrix33 &m, int p) {
	Matrix33 tmp(m);
	
	if(p == 0)
		return Matrix33();
	
	for(int i = 1; i < p; i++)
		tmp = multiply(tmp, m);
	
	return tmp;
}

inline Matrix33 rotMat3(Vector3 v, double angle){
	Matrix33 m;
	double ca = cos(angle);
	double sa = sin(angle);
	double one_ca = 1.0f - ca;

	m.elements[0][0] = (v.a * v.a) * one_ca + ca;
	m.elements[0][1] = (v.a * v.b) * one_ca - (v.c * sa);
	m.elements[0][2] = (v.a * v.c) * one_ca + (v.b * sa);

	m.elements[1][0] = (v.b * v.a) * one_ca + (v.c * sa);
	m.elements[1][1] = (v.b * v.b) * one_ca + ca;
	m.elements[1][2] = (v.b * v.c) * one_ca - (v.a * sa);

	m.elements[2][0] = (v.c * v.a) * one_ca - (v.b * sa);
	m.elements[2][1] = (v.c * v.b) * one_ca + (v.a * sa);
	m.elements[2][2] = (v.c * v.c) * one_ca + ca;

	return m;
}

inline double l_frobenius(Matrix33 &m){
	double sum = 0;
	sum += m.elements[0][0]*m.elements[0][0];
	sum += m.elements[0][1]*m.elements[0][1];
	sum += m.elements[0][2]*m.elements[0][2];
	sum += m.elements[1][0]*m.elements[1][0];
	sum += m.elements[1][1]*m.elements[1][1];
	sum += m.elements[1][2]*m.elements[1][2];
	sum += m.elements[2][0]*m.elements[2][0];
	sum += m.elements[2][1]*m.elements[2][1];
	sum += m.elements[2][2]*m.elements[2][2];
	return sqrt(sum);
}

inline double l_infinity(Matrix33 &m){
	double maxElem = 0;
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			double elem = fabs(m.elements[i][j]);
			if(elem > maxElem)
				maxElem = elem;
		}
	}

	return maxElem;
}

inline Matrix33 symmetricAdjoint(Matrix33 &M){
	double &a = M.elements[0][0];
	double &b = M.elements[0][1];
	double &d = M.elements[0][2];
	double &c = M.elements[1][1];
	double &e = M.elements[1][2];
	double &f = M.elements[2][2];

	double cf_e2 = c*f - e*e;
	double af_d2 = a*f - d*d;
	double ac_b2 = a*c - b*b;

	double bf_de = b*f - d*e;
	double be_dc = b*e - d*c;
	double ae_db = a*e - d*b;

	return Matrix33( cf_e2, -bf_de,  be_dc,
				    -bf_de,  af_d2, -ae_db,
				     be_dc, -ae_db,  ac_b2);
}

#endif
