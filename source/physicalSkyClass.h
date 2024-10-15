/*
    aaPhysicalSky Copyright(c) Amaan Akram.
	An Arnold shader for rendering sky colour and aerial perspective,
	based on "A Practical Analytic Model for Daylight" by Preetham et. al.
	Original paper: www.cs.utah.edu/~shirley/papers/sunsky/sunsky.pdf

    aaPhysicalSky is free software; you can redistribute it and/or modify
    it as long as you reproduce the the original author(s) names.

	Parts of this code are based on actual code from Brian Smits  
	(bes@phoenix.cs.utah.edu), co-author of the paper.

    aaPhysicalSky is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 */

#include <omp.h>
#include <string>
#include <math.h>
using namespace std;

#define USE_SSE2 1
#include "constants.h"
#include "fexp/t2exp.c"
#include "spectralData.h"
#include "functionLib.h"
#include "spectrum.h"
#include "spectrumGlobals.h"
#include "physicalSunClass.h"

class physicalSky
{
public:
	AtMutex my_mutex;
	bool is_init;
	
	Real T;				// turbitidy
	bool cieOvercast;	// sky model
	AtRGBA ground_color;	

	// angles in radians
	Real thetaSun;		// sun angle from zenith (radians)
	Real phiSun;		// sun angle from azimuth (radians) -- not needed for sky, and should be removed from here
	
	physicalSun pSun;	//physical Sun class object.

	//global colour multipliers
	Real multiplier;
	Real exposure;
	Real saturation;
	Real redBlueShift;
	Real candelaConversion; // Cd/m^2 to pixel intensity scaling factor

	//sun direction vector
	AtVector v_sunDir;

	// clear sky and overcast sky variables
	Real zenith_Y, zenith_x, zenith_y;
	Real perez_Y[5], perez_x[5], perez_y[5];

	//horizon data
	Real horizon_height;
	Real horizon_blur;

	//sun data
	Real sun_opacity;
	Real sun_intensity;
	Real sun_disk_scale;

	//multipliers for controlling sky luminance of different parts of the sky
	Real multA;
	Real multB;
	Real multC;
	Real multD;
	Real multE;

	// Aerial Perspective variables
	spectrum ST0_1[nTheta + 1][nPhi + 1];
	spectrum ST0_2[nTheta + 1][nPhi + 1];
	spectrum specOne;
	spectrum beta_m, beta_p, beta_m_angPrefix,  beta_p_angPrefix;
	Real fogTransparency, maxVisibilityDistance, atmosDepth;
	Real skyInscatter, sunInscatter; // multipliers for scattered radiance into viewer direction from sky and sun

	/////////
	Real maxTurbidity;
	Real alpha_1;			// Mie
	Real alpha_2;			// Rayleigh
	Real sunSolidAngle;
	Real anisoCorrection;  	// Correction for molecular anisotropy.
	Real N;					// Number of molecules per unit volume.
	
	AtVector zenith;
	AtVector azimuth;
	Real overcastWhitePoint_D65_x;
	Real overcastWhitePoint_D65_y;;
	/////////

	physicalSky();
	~physicalSky();

	// -----------------Defining methods below-----------------

	// Main Input function. A user needs to call this function to give input to this class's object
	void skyInput(Real, AtVector&,Real, Real, Real, Real, Real, AtRGBA&, 
				  Real, Real, Real, Real, Real, Real, Real, Real, Real, 
				  Real, Real, Real, Real, Real, bool);

	// -----------------Utility functions follow -- these should never be called directly-----------------

	// One-time init functions follow
	void reinit(Real Turbidity, AtVector &sunDir);
	void initSpectrums();
	void initDistributionCoEfficients();
	void initZenith();
	void sunDiscClamp(); // clamps the solar disc's intensity in the sky to make it easier to sample
	inline void computeAngles(Real& phi, Real& theta, const AtVector& v_dir);
	Real perez(int64_t x, Real gamma, Real theta, Real cosTheta, Real cos2gamma);

	// Aerial Perspective functions
	void ap_createConstants();
	void ap_InitS0();
	void ap_CalculateS0(Real thetav, Real phiv, spectrum &S0_1, spectrum &S0_2, AtVector& v_dir);
	spectrum ap_InscatteredRadiance(const Real h0, const Real s, AtVector& v_dir);
	spectrum ap_ExtinctionFactor(const Real h0, const Real s, AtVector& v_dir);
	spectrum ap_GetNeta(Real theta) ;
	void GetS0fromTable(Real theta, Real phi, spectrum &S0_1, spectrum &S0_2);
	
	// -----------------Per-sample functions follow-----------------

	// Adds horizon/ground 
	void addGround(AtRGB &skyCol, AtRGB &, const AtVector &v_camDir);
	// Adds sun in the sky
	void addSun(AtRGB &, AtRGB &, AtVector &gamma);
	// Main output function. This should be called to get colour of the sky per-sample of an environment shader
	// Output in candels / m^2
	spectrum getSkySpectralRadiance(AtVector &v_dir);
};

physicalSky::physicalSky()
{
	T = thetaSun = phiSun = multiplier = exposure = fogTransparency = redBlueShift = 0.0;
	horizon_height =  horizon_blur = sun_opacity = sun_intensity = sun_disk_scale = 0.0;
	multA = multB = multC = multD = multE = 0.0;
	ground_color.r = ground_color.g = ground_color.b = 0.0;
	v_sunDir.x = v_sunDir.y = v_sunDir.z = 0.0;
	maxVisibilityDistance = 500.0;
	atmosDepth = 32000.0;
	skyInscatter = sunInscatter = 1.0;
	cieOvercast = false;
	specOne = spectrum(1.0);

	maxTurbidity = 2.0;
	alpha_1 = 0.83331810215394e-3; // Mie
	alpha_2 = 0.11360016092149e-3; // Rayleigh
	sunSolidAngle = 0.25 * aa_PI * 1.39 * 1.39 / (150.0 * 150.0);
	anisoCorrection = 1.06;  	// Correction for molecular anisotropy.
	N = 2.545e25; 			// Number of molecules per unit volume.
	zenith.x = zenith.y = azimuth.y = azimuth.z = 0.0;
	zenith.z = azimuth.x = 1.0;
	overcastWhitePoint_D65_x = 0.3127f;
	overcastWhitePoint_D65_y = 0.3291f;
	is_init = false;
}

physicalSky::~physicalSky()
{
}

void physicalSky::skyInput(	
				Real		_turbidity,
				AtVector&	_v_sunDir,
				Real		_multA,
				Real		_multB,
				Real		_multC,
				Real		_multD,
				Real		_multE,
				AtRGBA&		_ground_color, 
				Real		_multiplier,
				Real		_horizon_height,
				Real		_horizon_blur,
				Real		_sun_opacity,
				Real		_sun_intensity,
				Real		_sun_disk_scale,
				Real		_saturation,
				Real		_candelaConversion,
				Real		_maxVisibilityDistance,
				Real		_fogTransparency,
				Real		_atmosDepth,
				Real		_sunInscatter,
				Real		_skyInscatter,
				Real		_redBlueShift,
				bool		_cieOvercast)
{
	multiplier			= maximum(_multiplier, 0.00001);
	ground_color.r		= _ground_color.r * _candelaConversion * 1.0/multiplier;  // we don't want to apply candela conversion to ground
	ground_color.g		= _ground_color.g * _candelaConversion * 1.0/multiplier;  // we don't want to apply candela conversion to ground
	ground_color.b		= _ground_color.b * _candelaConversion * 1.0/multiplier;  // we don't want to apply candela conversion to ground
	ground_color.a		= _ground_color.a;
	horizon_height		= _horizon_height * 0.01 * -10.0; // scaled for better GUI control. 0 is horizon. 100 is zenith;
	horizon_blur		= _horizon_blur   * 0.01; // scaled for better GUI control. 0 is horizon. 100 is zenith;
	sun_opacity			= _sun_opacity;
	sun_intensity		= _sun_intensity;
	sun_disk_scale		= _sun_disk_scale;
	saturation			= _saturation;
	candelaConversion	= 1.0 / _candelaConversion;
	maxVisibilityDistance = _maxVisibilityDistance;
	fogTransparency		= _fogTransparency;
	atmosDepth			= _atmosDepth;
	
	clamp(fogTransparency, 0.0, 1.0);
	bool reInitSky = 0;

	if(v_sunDir != _v_sunDir || T != _turbidity || sunInscatter != _sunInscatter || skyInscatter != _skyInscatter
		|| redBlueShift != _redBlueShift || cieOvercast != _cieOvercast)
	{
		sunInscatter	= _sunInscatter;
		skyInscatter	= _skyInscatter;
		redBlueShift	= _redBlueShift;
		cieOvercast		= _cieOvercast;
		reInitSky	= 1;
	}
	
	// no longer using multC. Am setting multC = multD later on.
	if( multA != _multA  || multB != _multB || multD != _multD 	|| multE != _multE ) 
	{
		multA = _multA;
		multB = _multB;
		multD = _multD;
		multE = _multE;
		reInitSky = 1;
	}
	if(reInitSky)
	{
		pSun.sunInput(_turbidity, _v_sunDir);
		reinit(_turbidity,_v_sunDir);
	}
}

void physicalSky::reinit(Real turbidity, AtVector &sunDir)
{
	T = turbidity; 
	v_sunDir = sunDir;

	computeAngles(phiSun, thetaSun, v_sunDir);

	initSpectrums();
	initDistributionCoEfficients();
	initZenith();
	ap_createConstants();
	ap_InitS0();
	sunDiscClamp();
}

void physicalSky::initSpectrums()
{
	#pragma omp parallel for
	for (int i = 0; i <= nTheta; ++i) 
	{
		for (int j = 0; j <= nPhi; ++j)
		{
			ST0_1[i][j] = spectrum(0.0);
			ST0_2[i][j] = spectrum(0.0);
		}
    }
}

void physicalSky::sunDiscClamp()
{
	// This function scales down the sun's intensity while maintaining its chromaticity information
	const Real maxSunLuminosity = 3.0;

	// get solar radiance of sun position in sky
	// not good to mirror sky
	AtVector skyDir = v_sunDir;
	skyDir.z = fabs(v_sunDir.z); 
	spectrum skySunSpot = getSkySpectralRadiance(skyDir);  // bad hack

	// convert solar radiance to xyY
	xyY skySunSpotxyY = skySunSpot.toxyY(); 
	// get xyY of attenuated sunlight
	XYZ sunSpotXYZ = RGBtoXYZ(pSun.sunDiscCol);
	xyY sunSpotxyY = XYZtoxyY(sunSpotXYZ);
	// clamp sun disc intensity to be at most 5 times brighter than sky
	clamp(sunSpotxyY.Y, aa_EPSILON, skySunSpotxyY.Y * maxSunLuminosity); 
	// store back in computed sun color
	xyYtoRGB(sunSpotxyY, pSun.sunDiscCol);  
}

inline void physicalSky::computeAngles(Real& phi, Real& theta, const AtVector& v_dir)
{
	Real dirz = v_dir.z;
	clamp(dirz,-1.0, 1.0);
	theta = acos(dirz);

	if(fabs(v_dir.x) < aa_EPSILON && fabs(v_dir.y) < aa_EPSILON)
		phi = 0.0;
	else
		phi = atan2(v_dir.y, v_dir.x);
}

void physicalSky::initDistributionCoEfficients()
{
	perez_Y[0] =  0.17872 * T - 1.46303;
	perez_Y[1] = -0.35540 * T + 0.42749;
	perez_Y[2] = -0.02266 * T + 5.32505;
	perez_Y[3] =  0.12064 * T - 2.57705;
	perez_Y[4] = -0.06696 * T + 0.37027;

	perez_x[0] = -0.01925 * T - 0.25922;
	perez_x[1] = -0.06651 * T + 0.00081;
	perez_x[2] = -0.00041 * T + 0.21247;
	perez_x[3] = -0.06409 * T - 0.89887;
	perez_x[4] = -0.00325 * T + 0.04517;

	perez_y[0] = -0.01669 * T - 0.26078;
	perez_y[1] = -0.09495 * T + 0.00921;
	perez_y[2] = -0.00792 * T + 0.21023;
	perez_y[3] = -0.04405 * T - 1.65369;
	perez_y[4] = -0.01092 * T + 0.05291;

	// See Page 8, last paragraph of Preetham et. al for explanation of these
	// A is the distribution coefficient for the darkening or brightening of the horizon
	// B is the distribution coefficient for the luminance gradient near the horizon
	// C is the distribution coefficient for the relative intensitiy of the circumsolar region (sun brightness)
	// D is the distribution coefficient for the width of the circumsolar region (sun size)
	// E is the distribution coefficient for relative backscattering
	// -----------------
	// not all of them seem to work intuitively. They all seem to have ranges of values over which they work

	// this seems to make the sun glow intensity work better
	// replaced multC with multD
	clamp(multD, 0.0, 2.0);

	perez_Y[0] *= multA;
	perez_Y[1] *= multB;
	perez_Y[2] *= multD;
	perez_Y[3] *= multD;
	perez_Y[4] *= multE;
}

void physicalSky::initZenith()
{
	Real theta2 = thetaSun * thetaSun;
	Real theta3 = theta2 * thetaSun;
	Real T2 = T * T;
	Real chi = (4.0 / 9.0 - T / 120.0) * (aa_PI - 2.0 * thetaSun);
	zenith_Y = (4.0453 * T - 4.9710)   * tan(chi) - 0.2155 * T + 2.4192;

	if(zenith_Y < 0.0) 
		zenith_Y =  aa_EPSILON;
	zenith_Y *= 1000.0;   // conversion from kilo cd/m^2 to cd/m^2

	zenith_x =
	    ( + 0.00165 * theta3 - 0.00374 * theta2 + 0.00208 * thetaSun) * T2 +
	    ( -0.02902  * theta3 + 0.06377 * theta2 - 0.03202 * thetaSun + 0.00394) * T +
	    ( + 0.11693 * theta3 - 0.21196 * theta2 + 0.06052 * thetaSun + 0.25885);

	zenith_y =
	    ( + 0.00275 * theta3 - 0.00610 * theta2 + 0.00316 * thetaSun) * T2 +
	    ( -0.04214  * theta3 + 0.08970 * theta2 - 0.04153 * thetaSun + 0.00515) * T +
	    ( + 0.15346 * theta3 - 0.26756 * theta2 + 0.06669 * thetaSun + 0.26688);

	chromaShift(zenith_x, zenith_y, redBlueShift);
}

Real physicalSky::perez(int64_t x, Real gamma, Real theta, Real cosTheta, Real cos2gamma)
{
	Real A, B, C, D, E;
	if(x == 1)
	{
		A = perez_x[0]; B = perez_x[1]; C = perez_x[2]; D = perez_x[3]; E = perez_x[4];
	}
	else if(x == 2)
	{
		A = perez_y[0]; B = perez_y[1]; C = perez_y[2]; D = perez_y[3]; E = perez_y[4];
	}
	else
	{
		A = perez_Y[0]; B = perez_Y[1]; C = perez_Y[2]; D = perez_Y[3]; E = perez_Y[4];
	}
	
	Real f_thetaGamma	= (1.0 + A * EXP( B /cosTheta)) * (1.0 + C * EXP(D * gamma) + E * cos2gamma);
	Real f_zeroThetaSun = (1.0 + A * EXP(B)) * (1.0 + C * EXP(D * thetaSun) + E * cos2gamma);

	// sometimes causes NaNs if not checked for division by zero
	if(fabs(f_zeroThetaSun) < aa_EPSILON) 
		return aa_EPSILON;
	else
		return f_thetaGamma / f_zeroThetaSun;
}

spectrum physicalSky::getSkySpectralRadiance(AtVector &v_dir)
{
	Real phi, theta, gamma, x, y, Y;

	computeAngles(phi,theta,v_dir);
	gamma = AiV3Dot(v_dir,v_sunDir);
	clamp(gamma,-1.0, 1.0);
	gamma = acos(gamma);

	Real cosTheta = cos(theta);
	Real cos2gamma = cos(gamma) * cos(gamma);

	if(!cieOvercast)
	{
		x = zenith_x * perez(1, gamma, theta, cosTheta, cos2gamma);
		y = zenith_y * perez(2, gamma, theta, cosTheta, cos2gamma);
		Y = zenith_Y * perez(3, gamma, theta, cosTheta, cos2gamma);
	}
	else
	{
		Real Y_occ = zenith_Y * ( 1.0 + cos(theta) / 3.0);
		x = overcastWhitePoint_D65_x; 
		y = overcastWhitePoint_D65_y;
		Y = zenith_Y * ( 1.0 + cos(theta) / 3.0);
	}

	// convert x, y chromaticities to spectral curve
	spectrum result(x, y);
	XYZ _XYZ = result.toXYZ();
	result = result * ( Y / _XYZ.Y);

	return result; 
}

void physicalSky::addGround(AtRGB &skyCol, AtRGB &inScatterCol, const AtVector &v_camDir)
{
	if(ground_color.a == 0.0f)
		return;

	AtRGB skyOriginal = skyCol;
	AtRGB grdColor;
	grdColor.r = ground_color.r +  inScatterCol.r;
	grdColor.g = ground_color.g +  inScatterCol.g;
	grdColor.b = ground_color.b +  inScatterCol.b;
	Real gamma = AiV3Dot(v_camDir,v_sunDir);
	clamp(gamma,-1.0, 1.0);
	gamma = acos(gamma);

	//set up horizon vector
	AtVector v_h = v_camDir;
	v_h = AiV3Normalize(v_camDir);
	v_h.z = horizon_height * -0.1; // for better UI control
	v_h = AiV3Normalize(v_h);

	//set up two vectors to point at horizon upper and lower blur limits
	AtVector v_h_blur_lower = v_h;
	v_h_blur_lower.z -= horizon_blur;
	v_h_blur_lower = AiV3Normalize(v_h_blur_lower);
	
	AtVector v_h_blur_uppper = v_h;
	v_h_blur_uppper.z += horizon_blur;
	v_h_blur_uppper = AiV3Normalize(v_h_blur_uppper);

	//check if current camera vector is pointing
	//somewhere in the upper and lower horizon blur region
	Real camDotH	= AiV3Dot(v_camDir, v_h_blur_lower);
	Real blurDotH   = AiV3Dot(v_h_blur_lower,v_h_blur_uppper);
	if( camDotH > blurDotH)
	{
		// linearly interpolate between sky colour and ground color
		Real t = rescale(camDotH, blurDotH, 1.0, 1.0 , 0.0, CLAMP_TO_RANGE);
		skyCol = skyCol *t*t + grdColor * (1.0 - t*t);
	}

	//add ground
	//do cross product to check for flipping ground direction
	AtVector v_cross;
	v_h_blur_lower = AiV3Cross(v_cross, v_camDir);
	if(v_camDir.x > 0.0)
	{
		if(v_cross.y < 0.0 )
			skyCol = grdColor;
	}
	else
	{
		if(v_cross.y > 0.0 )
			skyCol = grdColor;
	}
	if(ground_color.a < 1.0)
		skyCol = AiLerp( ground_color.a, skyOriginal, skyCol);

}

void physicalSky::addSun(AtRGB &skyCol, AtRGB &inScatterCol,  AtVector &v_dir)
{
	if(v_dir.z < 0.0 || cieOvercast) // sample ray going below horizon, or if overcast
		return;

	// physically correct sun size in radians = 0.5
	const Real sunSize = DegsToRads(0.5) * sun_disk_scale; 
	Real gamma = AiV3Dot(v_dir,v_sunDir);
	clamp(gamma,-1.0, 1.0);
	gamma = acos(gamma);
	const Real halfGamma = 0.5 * gamma;

	if( halfGamma < sunSize)
	{
		// multiplier to darken edges of the sun
		Real limbDarkening = pow(rescale(halfGamma, 0.0, sunSize, 1.0, 0.6), 0.5);

		// lerp factor to blend between sky and sun
		Real blurSun = pow(rescale(halfGamma, 0.0, sunSize, 1.0, 0.0), 0.35);

		Real opacity = minimum(sun_opacity, blurSun);

		skyCol = AiLerp(opacity, skyCol, (pSun.sunDiscCol * sun_intensity * limbDarkening));
	}
}

void physicalSky::ap_createConstants()
{
	initNetaTable();

	const Real c = (0.06544204545455 * T - 0.06509886363636) * 1e-15;
	const Real V = 4.0;
	const Real V2 = V - 2.0;

	// move var declarations out from here, and make private to omp pragma
	Real lambdasi, Nlambda4;
	Real n2_1;
	Real K;
	const Real _3N = 3.0 * N;
	const Real bm_mul1 = 8.0 * aa_PI * aa_PI * aa_PI * anisoCorrection;
	const Real bp_mul1 = 0.434 * aa_PI * c * pow(0.01, V - 3.0);
	const Real bmp_mul1 = aa_TWOPI * aa_PI * anisoCorrection * 0.7629;
	const Real bpp_mul1 = 0.217 * c * pow(0.01, V - 3.0);

	#pragma omp parallel for private(lambdasi, Nlambda4, n2_1, K)
	for (int i = 0; i < wavelengths; ++i)
	{
		lambdasi = spectral_data[i].wavelength * 1e-6;  /* Converstion to SI units (meters) */
		Nlambda4 = 1.0 / (_3N * lambdasi * lambdasi * lambdasi * lambdasi);

		/* Rayleigh total scattering coefficient */
		n2_1 = n2_1Amplitudes[i]; 
		beta_m.data[i] = bm_mul1 * n2_1 * Nlambda4;

		/* Mie total scattering coefficient */
		K = spectral_data[i].K;
		beta_p.data[i] = bp_mul1 * pow(aa_TWOPI/lambdasi, V2) * K;

		/* Rayleigh Angular scattering coefficient */
		beta_m_angPrefix.data[i] = bmp_mul1 * n2_1 * Nlambda4;  

		/* Mie Angular scattering coefficient */
		beta_p_angPrefix.data[i] =  bpp_mul1 * pow(aa_TWOPI / lambdasi, V2);
	}
}

spectrum physicalSky::ap_ExtinctionFactor(const Real h0, const Real s, AtVector& v_dir)
{
	Real theta, phi;
	computeAngles(phi,theta,v_dir);
	Real cosTheta = fabs(v_dir.z); // fabs to prevent negative theta, and hence negative attenuation below horizon
    Real B_1 = alpha_1 * cosTheta;
    Real B_2 = alpha_2 * cosTheta;
    Real constTerm_1 = EXP(-alpha_1 * h0) * EvalFunc(B_1, s);
    Real constTerm_2 = EXP(-alpha_2 * h0) * EvalFunc(B_2, s);
  
	constTerm_1 = -constTerm_1;
	constTerm_2 = -constTerm_2;

	spectrum a((beta_p * constTerm_1));
	spectrum b((beta_m * constTerm_2));
	a.exponentSpectrum();
	b.exponentSpectrum();
	a *= b;
	return a;
}

spectrum physicalSky::ap_InscatteredRadiance(const Real h0, const Real s, AtVector &v_dir) 
{
	Real theta, phi;
	computeAngles(phi,theta,v_dir);
	Real cosTheta = v_dir.z; 

	// to skip simpler approximation to avoid thin line near horizon
	if(fabs(cosTheta) < 0.001) 
	{
		if(cosTheta > 0.0)
			cosTheta = 0.001;
		else
			cosTheta = -0.001;
	}

    Real aCosTheta_1 = alpha_1 * cosTheta;
    Real aCosTheta_2 = alpha_2 * cosTheta;
	Real inv_aCosTheta_1 = 1.0 / aCosTheta_1;
    Real inv_aCosTheta_2 = 1.0 / aCosTheta_2;
    spectrum I_1, I_2;
  
	GetS0fromTable(theta, phi, I_1, I_2);

    // Analytical approximation
    Real A,B,C,D,H1,H2,K;
    Real u_f1, u_i1,u_f2, u_i2, int_f, int_i, fs, fdashs, fdash0;
    Real a1,b1,a2,b2;
    Real den1, den2;
	Real hsc = h0 + s * cosTheta;
    b1 = u_f1 = EXP(-alpha_1 * hsc);
    H1 = a1 = u_i1 = EXP(-alpha_1 * h0);
    b2 = u_f2 = EXP(-alpha_2 * hsc);
    H2 = a2 = u_i2 = EXP(-alpha_2 * h0);

	Real a1b1 = a1 - b1;
	Real a2b2 = a2 - b2;
    den1 = a1b1 * a1b1 * a1b1;
    den2 = a2b2 * a2b2 * a2b2;

	Real uf1 = u_i1 - u_f1;
	Real uf2 = u_i2 - u_f2;
    
    for (int i = 0; i < wavelengths; ++i) 
	{
		// for integral 1
		K = beta_p.data[i] * (inv_aCosTheta_1);
		fdash0 = -beta_m.data[i] * H2;
		fs = EXP(-beta_m.data[i] * (inv_aCosTheta_2) * uf2);
		fdashs = -fs * beta_m.data[i] * u_f2;

		CalculateABCD(a1, b1, fs, fdash0, fdashs, den1, A, B, C, D);
		int_f = Helper1(A, B, C, D, H1, K, u_f1);
		int_i = Helper1(A, B, C, D, H1, K, u_i1);
		I_1.data[i] *= ((int_f - int_i) * (-inv_aCosTheta_1));

		// for integral 2
		K = beta_m.data[i] * (inv_aCosTheta_2);
		fdash0 = -beta_p.data[i] * H1;
		fs = EXP(-beta_p.data[i] * (inv_aCosTheta_1) * uf1);
		fdashs = -fs * beta_p.data[i] * u_f1;

		CalculateABCD(a2, b2, fs, fdash0, fdashs, den2, A, B, C, D);
		int_f = Helper1(A, B, C, D, H2, K, u_f2);
		int_i = Helper1(A, B, C, D, H2, K, u_i2);
		I_2.data[i] *= ((int_f - int_i) * (-inv_aCosTheta_2));
    }
	I_1 += I_2;
    return I_1;
}

void physicalSky::ap_InitS0()
{
	int i, j;
    Real theta, phi;
	theta = phi =  0.0;
	
	AtVector zenith_dir = zenith;
	AtVector sample_dir = zenith;
    const Real delTheta =  180.0  / Real(nTheta);	// 9 degrees
    const Real delPhi   =  360.0  / Real(nPhi );	// 18 degrees
  
	// generate and rotate camera view directions
    for (i = 0; i <= nTheta; ++i) 
	{
		for (j = 0; j <= nPhi; ++j)
		{
			computeAngles(phi,theta, sample_dir);
			ap_CalculateS0(theta, phi,  ST0_1[i][j], ST0_2[i][j], sample_dir);
			sample_dir = vec_rotate_z(sample_dir,delPhi);
		}
		zenith_dir = vec_rotate_y(zenith_dir,delTheta);
		sample_dir = zenith_dir;
    }
}

void physicalSky::ap_CalculateS0(Real thetav, Real phiv, spectrum &S0_1, spectrum &S0_2, AtVector &v_Dir)
{
	Real delTheta = RadsToDegs(aa_PI / 2.0 / nTheta);	//  9 degrees
    Real delPhi   = RadsToDegs(aa_PI / nPhi);			// 18 degrees
	const Real delThetaPhi = delTheta * delPhi;
	Real theta, phi, psi, mul1, mul2;
	AtVector zenith_dir = zenith;
	AtVector sample_dir = zenith;

	spectrum skyRad;
	spectrum skyAmb_1, skyAmb_2;
	spectrum beta_ang_1, beta_ang_2;
	spectrum neta;

    for (int i = 0; i < nTheta; ++i) 
	{
		for (int j = 0; j < nPhi; ++j)
		{
			skyRad = getSkySpectralRadiance(sample_dir) * ap_ExtinctionFactor(1.f, 32000.0, sample_dir);

			computeAngles(phi,theta,sample_dir);
			psi = AiV3Dot(v_Dir, sample_dir);
			neta = ap_GetNeta(psi);
			mul1 = (1.0 + 0.9324 * cos(psi) * cos(psi));
			mul2 = sin(theta) * delThetaPhi;

			for(int k = 0; k < wavelengths; ++k)
			{
				beta_ang_1.data[k] = beta_p_angPrefix.data[k] * neta.data[k];
				beta_ang_2.data[k] = beta_m_angPrefix.data[k] * mul1;
	  
				skyAmb_1.data[k] += skyRad.data[k] * beta_ang_1.data[k] * mul2;
				skyAmb_2.data[k] += skyRad.data[k] * beta_ang_2.data[k] * mul2;
			}

			sample_dir = vec_rotate_z(sample_dir, delPhi);
		}
		zenith_dir = vec_rotate_y(zenith_dir, delTheta);
		sample_dir = zenith_dir;

	}
    /* Sun's ambience term*/
	psi = AiV3Dot(v_Dir, v_sunDir);
    beta_ang_1 = beta_p_angPrefix * ap_GetNeta(psi);
    beta_ang_2 = beta_m_angPrefix * (1 + 0.9324 * cos(psi) * cos(psi));
  
    spectrum sunAmb_1 = pSun.sunlightAttenuated * beta_ang_1;
    spectrum sunAmb_2 = pSun.sunlightAttenuated * beta_ang_2;
    
	const Real sunScale = 100.0;			// arbitrary multiplier. units issue?
	const Real skyScale = 1.0 / 10000.0;	// arbitrary multiplier. units issue?
	skyAmb_1 *= (skyInscatter * skyScale);
	skyAmb_2 *= (skyInscatter * skyScale);
	sunAmb_1 *= (sunScale * sunInscatter);
	sunAmb_2 *= (sunScale * sunInscatter);
	sunAmb_1 += skyAmb_1;
	sunAmb_2 += skyAmb_2;
    S0_1 += sunAmb_1; 
    S0_2 += sunAmb_2; 
}

// get spectrum for a specific direction
spectrum physicalSky::ap_GetNeta(Real theta)
{
	spectrum result;
	clamp(theta, -1.0, 1.0);
	theta = acos(theta);
	theta = RadsToDegs(theta) * 10.0;
	Real u = theta - (int)theta;
	if((theta < 0.0)||(theta > 1801.0)) 
	{
		AiMsgWarning("Theta outside range for GetNeta");
		return  result;
	}
	if (theta > 1800.0) 
		theta = 1800.0;

	if ((int)theta == 1800)
		return  netaTable[1800];

	 result =  netaTable[(int)theta] * (1-u)  +  netaTable[(int)theta+1] * u;
	 return result;

}

void physicalSky::GetS0fromTable(Real theta, Real phi, spectrum &S0_1, spectrum &S0_2)
{
	Real eps = 1e-4;
	if (phi < 0.0) 
		phi += 2.0 * aa_PI; // convert phi from -pi..pi to 0..2pi

	theta	= theta * nTheta / aa_PI - eps;
	phi		= phi * nPhi / (2.0 * aa_PI) - eps;
	
	if (theta < aa_EPSILON) 
		theta = 0.0;
	if (phi < aa_EPSILON) 
		phi = 0.0;

	int i = (int) theta;
	Real u = theta - i;
	int j = (int)phi;
	Real v = phi - j;

	S0_1 =	ST0_1[i][j] * ((1-u) * (1-v))	+ 
			ST0_1[i+1][j] * (u * (1-v))		+ 
			ST0_1[i][j+1] * ((1-u) * v)		+ 
			ST0_1[i+1][j+1] * u * v;

	S0_2 =	ST0_2[i][j] * (1-u) * (1-v)		+ 
			ST0_2[i+1][j] * u * (1-v)		+
			ST0_2[i][j+1] * (1-u) * v		+ 
			ST0_2[i+1][j+1] * u * v ;
}