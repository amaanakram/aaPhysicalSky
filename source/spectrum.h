#ifndef SPECTRUM
#define SPECTRUM

#include "colour.h"

class alignas(BOUNDARY) spectrum
{
public:
	// spectral data array
	alignas(BOUNDARY) Real data[wavelengths];

	// constructors
	spectrum();
	spectrum(Real init);
	spectrum(Real x, Real y);
	spectrum(Real normalArray[]);
	spectrum(const spectrum& copy);

	// operators
	spectrum operator + (const spectrum&) const;
	spectrum operator - (const spectrum&) const;
	spectrum operator * (const spectrum&) const;
	spectrum operator / (const spectrum&) const;
	spectrum operator + (const Real) const;
	spectrum operator - (const Real) const;
	spectrum operator * (const Real) const;
	spectrum operator / (const Real) const;

	spectrum& operator+=(const spectrum& in);
	spectrum& operator-=(const spectrum& in);
	spectrum& operator*=(const spectrum& in);
	spectrum& operator/=(const spectrum& in);

	spectrum& operator+=(const Real in);
	spectrum& operator-=(const Real in);
	spectrum& operator*=(const Real in);
	spectrum& operator/=(const Real in);

	spectrum operator-(void); 

	// functions
	void lerp(Real t, const spectrum& in);
	f_inline void exponentSpectrum();
	f_inline XYZ toXYZ();
	f_inline xyY toxyY();
	f_inline AtRGB toRGB();
};

spectrum::spectrum()
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		a_i = _mm_setzero_pd();
		_mm_store_pd(&data[i], a_i);
	}
}

spectrum::spectrum(Real init)
{
	v2d a_i =  _mm_set1_pd(init);
	for(int i = 0; i < wavelengths; i += 2)
		_mm_store_pd(&data[i], a_i);
}

spectrum::spectrum(Real x, Real y)
{
	Real M1 = (-1.3515 - 1.7703  * x +  5.9114 * y) / (0.0241 + 0.2562 * x - 0.7341 * y);
    Real M2 = ( 0.03   - 31.4424 * x + 30.0717 * y) / (0.0241 + 0.2562 * x - 0.7341 * y);

	v2d m1_i =  _mm_set1_pd(M1);
	v2d m2_i =  _mm_set1_pd(M2);
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d s0_i = _mm_load_pd(&spectral_data[i].S0);
		v2d s1_i = _mm_load_pd(&spectral_data[i].S1);
		v2d s2_i = _mm_load_pd(&spectral_data[i].S2);

		v2d s1m1_i = _mm_mul_pd(s1_i, m1_i);
		v2d s2m2_i = _mm_mul_pd(s2_i, m2_i);

		v2d out_i = _mm_add_pd(s0_i, s1m1_i);
		out_i = _mm_add_pd(out_i, s2m2_i);
		_mm_store_pd(&data[i], out_i);
	}
}

spectrum::spectrum(Real normalArray[]) // no range checks!
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&normalArray[i]);
		_mm_store_pd(&data[i], a_i);
	}
}

spectrum::spectrum(const spectrum& copy) 
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&copy.data[i]);
		_mm_store_pd(&data[i], a_i);
	}
}

spectrum spectrum::operator+ (const spectrum& in) const
{
	spectrum result;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&in.data[i]);
		v2d b_i = _mm_load_pd(&data[i]);
 
		a_i = _mm_add_pd(a_i, b_i);
 
		_mm_store_pd(&result.data[i], a_i);
	}
	return result;

}

spectrum spectrum::operator- (const spectrum& in) const
{
	spectrum result;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&in.data[i]);
		v2d b_i = _mm_load_pd(&data[i]);
 
		a_i = _mm_sub_pd(a_i, b_i);
 
		_mm_store_pd(&result.data[i], a_i);
	}
	return result;
}

spectrum spectrum::operator* (const spectrum& in) const
{
	spectrum result;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&in.data[i]);
		v2d b_i = _mm_load_pd(&data[i]);
 
		a_i = _mm_mul_pd(a_i, b_i);
 
		_mm_store_pd(&result.data[i], a_i);
	}
	return result;
}

spectrum spectrum::operator/ (const spectrum& in) const
{
	spectrum result;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&in.data[i]);
		v2d b_i = _mm_load_pd(&data[i]);
 
		a_i = _mm_div_pd(a_i, b_i);
 
		_mm_store_pd(&result.data[i], a_i);
	}
	return result;
}

spectrum spectrum::operator+ (const Real in) const
{
	spectrum result;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_set1_pd(in);
		v2d b_i = _mm_load_pd(&data[i]);
 
		a_i = _mm_add_pd(a_i, b_i);
 
		_mm_store_pd(&result.data[i], a_i);
	}
	return result;
}

spectrum spectrum::operator- (const Real in)  const
{
	spectrum result;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_set1_pd(in);
		v2d b_i = _mm_load_pd(&data[i]);
 
		a_i = _mm_sub_pd(a_i, b_i);
 
		_mm_store_pd(&result.data[i], a_i);
	}
	return result;
}

spectrum spectrum::operator* (const Real in) const
{
	spectrum result;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_set1_pd(in);
		v2d b_i = _mm_load_pd(&data[i]);
 
		a_i = _mm_mul_pd(a_i, b_i);
 
		_mm_store_pd(&result.data[i], a_i);
	}
	return result;
}

spectrum spectrum::operator/ (const Real in) const
{
	spectrum result;
	const Real inv_in = 1.0 / in;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_set1_pd(inv_in);
		v2d b_i = _mm_load_pd(&data[i]);
 
		a_i = _mm_mul_pd(a_i, b_i);
 
		_mm_store_pd(&result.data[i], a_i);
	}
	return result;
}

spectrum& spectrum::operator+=(const spectrum& in)
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		v2d b_i = _mm_load_pd(&in.data[i]);
		a_i = _mm_add_pd(a_i, b_i);
		_mm_store_pd(&data[i], a_i);
	}
	return *this;
}
spectrum& spectrum::operator-=(const spectrum& in)
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		v2d b_i = _mm_load_pd(&in.data[i]);
		a_i = _mm_sub_pd(a_i, b_i);
		_mm_store_pd(&data[i], a_i);
	}
	return *this;
}
spectrum& spectrum::operator*=(const spectrum& in)
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		v2d b_i = _mm_load_pd(&in.data[i]);
		a_i = _mm_mul_pd(a_i, b_i);
		_mm_store_pd(&data[i], a_i);
	}
	return *this;
}
spectrum& spectrum::operator/=(const spectrum& in)
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		v2d b_i = _mm_load_pd(&in.data[i]);
		a_i = _mm_div_pd(a_i, b_i);
		_mm_store_pd(&data[i], a_i);
	}
	return *this;
}

spectrum& spectrum::operator+=(const Real in)
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		v2d b_i = _mm_set1_pd(in);
		a_i = _mm_add_pd(a_i, b_i);
		_mm_store_pd(&data[i], a_i);
	}
	return *this;
}
spectrum& spectrum::operator-=(const Real in)
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		v2d b_i = _mm_set1_pd(in);
		a_i = _mm_sub_pd(a_i, b_i);
		_mm_store_pd(&data[i], a_i);
	}
	return *this;
}
spectrum& spectrum::operator*=(const Real in)
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		v2d b_i = _mm_set1_pd(in);
		a_i = _mm_mul_pd(a_i, b_i);
		_mm_store_pd(&data[i], a_i);
	}
	return *this;
}
spectrum& spectrum::operator/=(const Real in)
{
	const Real inv = 1.0 / in;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		v2d b_i = _mm_set1_pd(inv);
		a_i = _mm_mul_pd(a_i, b_i);
		_mm_store_pd(&data[i], a_i);
	}
	return *this;
}

spectrum spectrum::operator-(void)
{
	spectrum result;
	for(int i = 0; i < wavelengths; i += 2)
	{
		v2d a_i = _mm_load_pd(&data[i]);
		v2d b_i = _mm_set1_pd(-1.0);
		a_i = _mm_mul_pd(a_i, b_i);
 
		_mm_store_pd(&result.data[i], a_i);
	}
	return result;
}

f_inline void spectrum::exponentSpectrum ()
{
	for(int i = 0; i < wavelengths; i += 2)
	{
		this->data[i] = EXP(this->data[i]);
		this->data[i+1] = EXP(this->data[i+1]);
	}
}

void spectrum::lerp(Real t, const spectrum& in)
{
	if(t <= 1.0)
	{
		v2d t_i  = _mm_set1_pd(t);
		v2d t1_i = _mm_set1_pd(1.0 - t);

		for(int i = 0; i < wavelengths; i += 2)
		{
			v2d a_i = _mm_load_pd(&data[i]);
			v2d b_i = _mm_load_pd(&in.data[i]);

			a_i = _mm_mul_pd(a_i, t_i);
			b_i = _mm_mul_pd(b_i, t1_i);
			a_i = _mm_add_pd(a_i, b_i);

			_mm_store_pd(&data[i], a_i);
		}
	}
	else
	{
		for(int i = 0; i < wavelengths; i += 2)
			_mm_store_pd(&data[i], _mm_load_pd(&in.data[i]));
	}
}

f_inline XYZ spectrum::toXYZ()
{
	XYZ _XYZ;

	for (int i = 0; i < wavelengths; i += 2) 
	{
		_XYZ.X += this->data[i] * cie_table[i][4] * wavelengthDelta; 
        _XYZ.Y += this->data[i] * cie_table[i][5] * wavelengthDelta; 
        _XYZ.Z += this->data[i] * cie_table[i][6] * wavelengthDelta;

		_XYZ.X += this->data[i+1] * cie_table[i+1][4] * wavelengthDelta; 
        _XYZ.Y += this->data[i+1] * cie_table[i+1][5] * wavelengthDelta; 
        _XYZ.Z += this->data[i+1] * cie_table[i+1][6] * wavelengthDelta;
    }
	return _XYZ;
}

f_inline xyY spectrum::toxyY()
{
	XYZ _XYZ = toXYZ();
	xyY _xyY = XYZtoxyY(_XYZ);
	return _xyY;
}

f_inline AtRGB spectrum::toRGB()
{
	XYZ _XYZ = toXYZ();
	AtRGB result = XYZtoRGB(_XYZ);
	return result;
}

#endif // SPECTRUM_H