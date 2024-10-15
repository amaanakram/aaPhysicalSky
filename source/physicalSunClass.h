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

#include "constants.h"

class physicalSun
{
public:
	Real thetaSun; // sun angle from zenith (radians)
	Real phiSun;   // sun angle from azimuth (radians)
	Real sunSolidAngle;
	AtVector v_sunDir;			 // sun direction vector
	Real T;						 // Turbidity for sunlight attenuation
	AtRGB sunLightCol;			 // colour for light shader
	AtRGB sunDiscCol;			 // colour for drawing sun disc
	Real candelaConversion;		 // not used
	spectrum sunlightAttenuated; // final sunlight spectrum on earth's surface

	physicalSun()
	{
		thetaSun = phiSun = 0.0;
		T = candelaConversion = 1.0;
		v_sunDir.x = v_sunDir.y = v_sunDir.z = 0.0;
		sunLightCol.r = sunLightCol.g = sunLightCol.b = 0.0;
		sunSolidAngle = 0.25 * aa_PI * 1.39 * 1.39 / (150.0 * 150.0);
	}
	~physicalSun()
	{
	}

	// utility function:theta and phi calculation function
	inline void computeAngles(const AtVector &v_dir);

	// The following two functions should be called from wherever this class is needed
	// Examples: an environment shader, or a physically correct light shader

	// Main input function
	void sunInput(Real _turbidity, AtVector &_v_sunDir);

	// Main output function - computes and returns the final sun attentuated spectrum
	void computeAttenuatedSunlight();
	spectrum getAttenuatedSunlight();
};

spectrum physicalSun::getAttenuatedSunlight()
{
	return sunlightAttenuated;
}

void physicalSun::sunInput(Real _turbidity, AtVector &_v_sunDir)
{
	if (T != _turbidity || v_sunDir != _v_sunDir)
	{
		T = _turbidity;

		v_sunDir = _v_sunDir;
		computeAngles(v_sunDir);
		computeAttenuatedSunlight();
		// convert to XYZ model
		XYZ _XYZ = sunlightAttenuated.toXYZ();

		// return result in cd/m^2
		// see en.wikipedia.org/wiki/Candela for the 683 scale factor (CIE Ybar has been divided by 683 to normalize)
		xyY _xyY = XYZtoxyY(_XYZ);
		_xyY.Y *= 683.002;
		xyYtoXYZ(_xyY, _XYZ);
		sunLightCol = XYZtoRGB(_XYZ);

		_xyY.Y /= sunSolidAngle;
		xyYtoXYZ(_xyY, _XYZ);
		sunDiscCol = XYZtoRGB(_XYZ);
	}
}

void physicalSun::computeAngles(const AtVector &v_dir)
{
	Real dirz = v_dir.z;
	clamp(dirz, -1.0, 1.0);
	thetaSun = acos(dirz);

	if (fabs(v_dir.x) < aa_EPSILON && fabs(v_dir.y) < aa_EPSILON)
		phiSun = 0.0;
	else
		phiSun = atan2(v_dir.y, v_dir.x);
}

void physicalSun::computeAttenuatedSunlight()
{
	// Working variables
	Real tauR, tauA, tauO, tauWA;
	Real beta = 0.04608365822050 * T - 0.04586025928522; // turbidity coefficient
	const Real alpha = 1.3;								 // wavelength exponent
	const Real lOzone = 0.35;							 // ozone at NTP (cm); NTP = Normal Temperature and Pressure
	const Real w = 2.0;									 // precipitable water vapor (cm)

	// prevent theta greater than 90 degrees to prevent optical mass going to infinity
	clamp(thetaSun, 0.0, aa_PIBYTWO);

	// Relative Optical Mass
	const Real m = (1.0 / (cos(thetaSun) + 0.15 * pow(Real(93.885 - thetaSun * aa_180BYPI), -1.253)));

	const Real tauRmul = -m * 0.008735;
	const Real tauAmul = -m * beta;
	const Real tauOmul = -m * lOzone;

	// We compute attenuation/absorption coefficients per wavelength band
	// Sun's spectral radiance is attenuated/scaled with these absoprtion coefficients per wavelength band

#pragma omp parallel for
	for (int i = 0; i < wavelengths; ++i)
	{
		// Rayleigh Scattering
		tauR = EXP(tauRmul * pow(spectral_data[i].wavelength, -4.08));

		// Aerosal (water + dust) attenuation
		tauA = EXP(tauAmul * pow(spectral_data[i].wavelength, -alpha));

		// Attenuation due to ozone absorption
		if (spectral_data[i].k_o > 0.0)
			tauO = EXP(tauOmul * spectral_data[i].k_o);
		else
			tauO = 1.0;

		// Attenuation due to water vapor absorbtion
		if (spectral_data[i].k_wa > 0.0)
			tauWA = EXP(-0.2385 * spectral_data[i].k_wa * w * m / pow(1.0 + 20.07 * spectral_data[i].k_wa * w * m, 0.45));
		else
			tauWA = 1.0;

		// populate spectral array
		sunlightAttenuated.data[i] = spectral_data[i].sun_spectral_radiance * tauR * tauA * tauO * tauWA;
		sunlightAttenuated.data[i] *= 10000.0;		 // convert from per-cm^2 to per-m^2
		sunlightAttenuated.data[i] /= 1000.0;		 // convert from per-micrometer to per-nanometer;
		sunlightAttenuated.data[i] *= sunSolidAngle; // convert from per steradian to per sun's solid angle
	}
}