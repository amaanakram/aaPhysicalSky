#ifndef COLOUR
#define COLOUR

// CIE spectral sensitivy curves
// Wavelenght in nm, 380 to 780nm
// using 2 and 10 degree observer CIE functions from
// www.cie.co.at/publ/abst/datatables15_2004/x2.txt
// www.cie.co.at/publ/abst/datatables15_2004/y2.txt
// www.cie.co.at/publ/abst/datatables15_2004/z2.txt
// www.cie.co.at/publ/abst/datatables15_2004/x10.txt
// www.cie.co.at/publ/abst/datatables15_2004/y10.txt
// www.cie.co.at/publ/abst/datatables15_2004/z10.txt


static const Real cie_table[41][7] =
{//lamda	x2bar		y2bar		z2bar		x10bar		y10bar		z10bar
{380,		0.001368,	0.000039,	0.006450,	0.000160,	0.000017,	0.000705},
{390,		0.004243,	0.000120,	0.020050,	0.002362,	0.000253,	0.010482},
{400,		0.014310,	0.000396,	0.067850,	0.019110,	0.002004,	0.086011},
{410,		0.043510,	0.001210,	0.207400,	0.084736,	0.008756,	0.389366},
{420,		0.134380,	0.004000,	0.645600,	0.204492,	0.021391,	0.972542},
{430,		0.283900,	0.011600,	1.385600,	0.314679,	0.038676,	1.553480},
{440,		0.348280,	0.023000,	1.747060,	0.383734,	0.062077,	1.967280},
{450,		0.336200,	0.038000,	1.772110,	0.370702,	0.089456,	1.994800},
{460,		0.290800,	0.060000,	1.669200,	0.302273,	0.128201,	1.745370},
{470,		0.195360,	0.090980,	1.287640,	0.195618,	0.185190,	1.317560},
{480,		0.095640,	0.139020,	0.812950,	0.080507,	0.253589,	0.772125},
{490,		0.032010,	0.208020,	0.465180,	0.016172,	0.339133,	0.415254},
{500,		0.004900,	0.323000,	0.272000,	0.003816,	0.460777,	0.218502},
{510,		0.009300,	0.503000,	0.158200,	0.037465,	0.606741,	0.112044},
{520,		0.063270,	0.710000,	0.078250,	0.117749,	0.761757,	0.060709},
{530,		0.165500,	0.862000,	0.042160,	0.236491,	0.875211,	0.030451},
{540,		0.290400,	0.954000,	0.020300,	0.376772,	0.961988,	0.013676},
{550,		0.433450,	0.994950,	0.008750,	0.529826,	0.991761,	0.003988},
{560,		0.594500,	0.995000,	0.003900,	0.705224,	0.997340,	0.000000},
{570,		0.762100,	0.952000,	0.002100,	0.878655,	0.955552,	0.000000},
{580,		0.916300,	0.870000,	0.001650,	1.014160,	0.868934,	0.000000},
{590,		1.026300,	0.757000,	0.001100,	1.118520,	0.777405,	0.000000},
{600,		1.062200,	0.631000,	0.000800,	1.123990,	0.658341,	0.000000},
{610,		1.002600,	0.503000,	0.000340,	1.030480,	0.527963,	0.000000},
{620,		0.854450,	0.381000,	0.000190,	0.856297,	0.398057,	0.000000},
{630,		0.642400,	0.265000,	0.000050,	0.647467,	0.283493,	0.000000},
{640,		0.447900,	0.175000,	0.000020,	0.431567,	0.179828,	0.000000},
{650,		0.283500,	0.107000,	0.000000,	0.268329,	0.107633,	0.000000},
{660,		0.164900,	0.061000,	0.000000,	0.152568,	0.060281,	0.000000},
{670,		0.087400,	0.032000,	0.000000,	0.081261,	0.031800,	0.000000},
{680,		0.046770,	0.017000,	0.000000,	0.040851,	0.015905,	0.000000},
{690,		0.022700,	0.008210,	0.000000,	0.019941,	0.007749,	0.000000},
{700,		0.011359,	0.004102,	0.000000,	0.009577,	0.003718,	0.000000},
{710,		0.005790,	0.002091,	0.000000,	0.004553,	0.001768,	0.000000},
{720,		0.002899,	0.001047,	0.000000,	0.002175,	0.000846,	0.000000},
{730,		0.001440,	0.000520,	0.000000,	0.001045,	0.000407,	0.000000},
{740,		0.000690,	0.000249,	0.000000,	0.000508,	0.000199,	0.000000},
{750,		0.000332,	0.000120,	0.000000,	0.000251,	0.000098,	0.000000},
{760,		0.000166,	0.000060,	0.000000,	0.000126,	0.000050,	0.000000},
{770,		0.000083,	0.000030,	0.000000,	0.000065,	0.000025,	0.000000},
{780,		0.000042,	0.000015,	0.000000,	0.000033,	0.000013,	0.000000}
};

struct alignas(BOUNDARY) XYZ
{
	Real X,Y,Z;

	XYZ() : X(0.0), Y(0.0), Z(0.0) {}
    XYZ(Real _X, Real _Y, Real _Z) : X(_X), Y(_Y), Z(_Z) {}
};

struct alignas(BOUNDARY) xyY
{
	Real x,y,Y;

	xyY() : x(0.0), y(0.0), Y(0.0) {}
    xyY(Real _x, Real _y, Real _Y) : x(_x), y(_y), Y(_Y) {}
};

struct alignas(BOUNDARY) xyz
{
	Real x,y,z;

    xyz() : x(0.0), y(0.0), z(0.0) {}
    xyz(Real _x, Real _y, Real _z) : x(_x), y(_y), z(_z) {}
};

void clampRGB(AtRGB& dest)
{
	if(dest.r < aa_EPSILON)
		dest.r = 0.0;
	if(dest.g < aa_EPSILON)
		dest.g = 0.0;
	if(dest.b < aa_EPSILON)
		dest.b = 0.0;
}

f_inline AtRGB XYZtoRGB(const XYZ &_XYZ)
{
	// CIE  XYZ into RGB using sRGB primaries and D65 white point
	// http://www.poynton.com/notes/colour_and_gamma/ColorFAQ.html
	// see "How do I transform between CIE  XYZ and a particular set of RGB primaries?"
	AtRGB dest;
	dest.r =  3.2404542f * _XYZ.X - 1.5371385f * _XYZ.Y - 0.4985314f * _XYZ.Z;
	dest.g = -0.9692660f * _XYZ.X + 1.8760108f * _XYZ.Y + 0.0415560f * _XYZ.Z;
	dest.b =  0.0556434f * _XYZ.X - 0.2040259f * _XYZ.Y + 1.0572252f * _XYZ.Z;	

	clampRGB(dest);
	return dest;
}

XYZ RGBtoXYZ(AtRGB color)
{
	XYZ result;
	result.X = 0.4124564 * color.r +  0.3575761 * color.g + 0.1804375 * color.b;
	result.Y = 0.2126729 * color.r +  0.7151522 * color.g + 0.0721750 * color.b;
	result.Z = 0.0193339 * color.r +  0.1191920 * color.g + 0.9503041 * color.b;
	return result;
}

void xyYtoXYZ( const xyY &src, XYZ &dest )
{
	if(src.y > aa_EPSILON)
	{
		dest.X = src.x * (src.Y / src.y);
		dest.Y = src.Y;
		dest.Z = (1.0 - src.x - src.y) * (src.Y / src.y);
	}
	else
		dest.X = dest.Y = dest.Z = 0.0;
}

xyY XYZtoxyY( const XYZ &src)
{
	xyY _xyY;
	Real XYZ = (src.X + src.Y + src.Z);
	if (XYZ > aa_EPSILON)
	{
		_xyY.x = src.X / XYZ;
		_xyY.y = src.Y / XYZ;
		_xyY.Y = src.Y;
	}
	return _xyY;
}

void xyYtoRGB( const xyY &src, AtRGB &dest)
{
	XYZ _XYZ;
	xyYtoXYZ(src, _XYZ);
	dest = XYZtoRGB(_XYZ);
}

xyz XYZtoxyz(XYZ& _XYZ)
{
	xyz _xyz;
	_xyz.x = _XYZ.X / (_XYZ.X + _XYZ.Y + _XYZ.Z);
	_xyz.y = _XYZ.Y / (_XYZ.X + _XYZ.Y + _XYZ.Z);
	_xyz.z = 1.0 - _xyz.x - _xyz.y;
	return _xyz;
}

void RGBtoHSV(const Real r,const Real g,const Real b, Real *h, Real *s, Real *v )
{
	Real min, max, delta;
	min = AiMin( r, g, b );
	max = AiMax( r, g, b );
	*v = max;
	delta = max - min;
	if( max != 0 )
		*s = delta / max;
	else 
	{
		*s = 0;
		*h = -1;
		return;
	}
	if( r == max )
		*h = ( g - b ) / delta;		// between yellow & magenta
	else if( g == max )
		*h = 2 + ( b - r ) / delta;	// between cyan & yellow
	else
		*h = 4 + ( r - g ) / delta;	// between magenta & cyan
	*h *= 60; // degrees
	if( *h < 0 )
		*h += 360;
}

void HSVtoRGB( Real *r, Real *g, Real *b, Real h, Real s, Real v )
{
	int i;
	Real f, p, q, t;
	if( s == 0 ) 
	{
		*r = *g = *b = v;
		return;
	}
	h /= 60; // sector 0 to 5
	i = floor( h );
	f = h - i; // factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );
	switch( i ) 
	{
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;
		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;
		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;
		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;
		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;
		default:		// case 5:
			*r = v;
			*g = p;
			*b = q;
			break;
	}
}

void chromaShift(Real& zen_x, Real& zen_y, Real t)
{
	// shit zenith x, y chromaticity
	// custom x,y primaries
	const Real red_x = 0.52;
	const Real red_y = 0.41;

	const Real blue_x = 0.16;
	const Real blue_y = 0.1;

	// going the full range produces unwanted results
	//t = rescale(t,-1.0, 1.0, -0.4, 0.4, CLAMP_TO_RANGE);

	if(t < 0.0) // red tint
	{
		zen_x = lerp(fabs(t), zen_x, red_x);
		zen_y = lerp(fabs(t), zen_y, red_y);
	}
	else	// blue tint
	{
		zen_x = lerp(t, zen_x, blue_x);
		zen_y = lerp(t, zen_y, blue_y);
	}
}

#endif // COLOUR_H

