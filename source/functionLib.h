#ifndef FUNCTIONLIB
#define FUNCTIONLIB

#include "constants.h"

f_inline Real DegsToRads(Real degrees) // Degrees to radians conversion...
{
	return (degrees * aa_PIBY180);
}

f_inline Real RadsToDegs(Real rads) // Radians to degrees conversion...
{
	return (rads * aa_180BYPI);
}

inline Real lerp(Real t, Real a, Real b)
{
	// t = 0 retunes a
	// t = 1 returns b
	return a * (1.0 - t) + b * t;
}

AtVector AnglesToVector(Real theta, Real phi)
{
	AtVector v;

	theta = DegsToRads(theta);
	phi = DegsToRads(180.0 - phi);

	v.x = cos(phi) * cos(theta);
	v.y = sin(phi) * cos(theta);
	v.z = sin(theta);

	return v;
}

inline void clamp(Real &x, Real a, Real b)
{
	if (x < a)
		x = a;
	if (x > b)
		x = b;
}

Real rescale(const Real &value, const Real &oldMin, const Real &oldMax, const Real &newMin, const Real &newMax, bool doClamp = 0)
{
	const Real oldDistance = oldMax - oldMin;
	const Real newDistance = newMax - newMin;
	const Real distance = (value - oldMin) / oldDistance;
	Real newValue = newMin + (distance * newDistance);
	if (doClamp)
	{
		if (newMin < newMax)
			clamp(newValue, newMin, newMax);
		else
			clamp(newValue, newMax, newMin);
	}
	return newValue;
}

bool almostEqual(Real A, Real B, Real maxRelativeError)
{
	if (A == B)
		return true;
	Real relativeError = fabs((A - B) / B);
	if (relativeError <= maxRelativeError)
		return true;
	return false;
}

template <class T>
const T &maximum(const T &a, const T &b)
{
	return (b < a) ? a : b;
}

template <class T>
const T &minimum(const T &a, const T &b)
{
	return (b > a) ? a : b;
}

f_inline Real EvalFunc(Real alpha_CosTheta, Real x) // from Smits
{
	if (fabs(alpha_CosTheta * x) < 0.01)
		return x;
	return (1.0 - EXP(-alpha_CosTheta * x)) / alpha_CosTheta;
}

f_inline Real Helper1(Real A, Real B, Real C, Real D, Real H, Real K, Real u)
{
	Real u2 = u * u;
	Real K2 = 1.0 / (K * K);
	Real inv_K = 1.0 / K;

	Real t = EXP(-K * (H - u));
	return (t * inv_K) * ((A * u * u2 + B * u2 + C * u + D) -
						  (3.0 * A * u2 + (B + B) * u + C) * inv_K +
						  (6.0 * A * u + B + B) * (K2) -
						  (6.0 * A) * (inv_K * K2));
}

inline void CalculateABCD(Real a, Real b, Real c, Real d, Real e,
						  Real den, Real &A, Real &B, Real &C, Real &D)
{
	Real inv_den = 1.0 / den;
	Real a3 = a * a * a;
	Real a2 = a * a;
	Real b3 = b * b * b;
	Real b2 = b * b;
	Real bd = b * d;
	Real ab = a * b;
	Real ad = a * d;
	Real ae = a * e;
	Real _3a = 3.0 * a;
	Real _3b = 3.0 * b;
	Real _6ab = 6.0 * ab;

	A = (-bd - 2.0 + c + c + ae - b * e + ad) * inv_den;

	B = -((a2 + a2) * e + a2 * d - _3a - ab * e + _3a * c + a * bd - (b2 + b2) * d - _3b - b2 * e + _3b * c) * inv_den;

	C = (-b3 * d - (b2 + b2) * ae - b2 * ad + ae * ab + (a2 + a2) * bd - _6ab + _6ab * c + a3 * e) * inv_den;

	D = -(b3 - b3 * ad - b2 * a2 * e + b2 * a2 * d - _3a * b2 + b * e * a3 - c * a3 + _3b * c * a2) * inv_den;
}

AtVector vec_rotate_x(AtVector &v_dir, Real theta)
{
	AtVector result;
	theta = DegsToRads(theta);
	result.x = v_dir.x;
	result.y = v_dir.y * cos(theta) - v_dir.z * sin(theta);
	result.z = v_dir.y * sin(theta) + v_dir.z * cos(theta);
	return result;
}

AtVector vec_rotate_y(AtVector &v_dir, Real theta)
{
	AtVector result;
	theta = DegsToRads(theta);
	result.x = v_dir.x * cos(theta) + v_dir.z * sin(theta);
	result.y = v_dir.y;
	result.z = v_dir.z * cos(theta) - v_dir.x * sin(theta);
	return result;
}

AtVector vec_rotate_z(AtVector &v_dir, Real theta)
{
	AtVector result;
	theta = DegsToRads(theta);
	result.x = v_dir.x * cos(theta) - v_dir.y * sin(theta);
	result.y = v_dir.x * sin(theta) + v_dir.y * cos(theta);
	result.z = v_dir.z;
	return result;
}

#endif //  FUNCTIONLIB_H