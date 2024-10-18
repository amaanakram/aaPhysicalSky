#include <string.h>
#include "ai.h"
#include <algorithm>
#include <cstddef> // for std::size_t if needed elsewhere
#include <iostream>
#include <cmath>

#include "physicalSkyClass.h"

static float far_clip; // camera far clipping distance

AI_SHADER_NODE_EXPORT_METHODS(aaPhysicalSkyMethods);

enum
{
	p_on,
	p_y_is_up,
	p_turbidity,
	p_sun_dir_x,
	p_sun_dir_y,
	p_sun_dir_z,
	p_multiplier,
	p_ground_color,
	p_horizon_height,
	p_horizon_blur,
	p_sun_opacity,
	p_sun_disk_scale,
	p_sun_intensity,
	p_mult_a,
	p_mult_b,
	p_mult_c,
	p_mult_d,
	p_mult_e,
	p_saturation,
	p_ray_reflected,
	p_ray_refracted,
	p_ray_diffuse,
	p_ray_glossy,
	p_ray_camera,
	p_candela_conversion,
	p_max_visibility_distance,
	p_fog_transparency,
	p_atmos_depth,
	p_sun_inscatter,
	p_sky_inscatter,
	p_red_blue_shift,
	p_tonemap,
	p_inscatter_fade,
	p_attenuation_fade,
	p_cie_overcast,
	p_log_sun_intensity,
	p_camera_height
};

node_parameters
{
	AiParameterBool("on", 1);
	AiParameterBool("y_is_up", 1);
	AiParameterFlt("turbidity", 5.0f);
	AiParameterFlt("sun_dir_x", 0.0f);
	AiParameterFlt("sun_dir_y", 0.0f);
	AiParameterFlt("sun_dir_z", 0.0f);
	AiParameterFlt("multiplier", 1.0f);
	AiParameterRGBA("ground_color", 0.2f, 0.2f, 0.2f, 1.0f);
	AiParameterFlt("horizon_height", 0.0f);
	AiParameterFlt("horizon_blur", 0.1f);
	AiParameterFlt("sun_opacity", 1.0f);
	AiParameterFlt("sun_disk_scale", 1.0f);
	AiParameterFlt("sun_intensity", 10.0f);
	AiParameterFlt("mult_a", 1.0f);
	AiParameterFlt("mult_b", 1.0f);
	AiParameterFlt("mult_c", 1.0f);
	AiParameterFlt("mult_d", 1.0f);
	AiParameterFlt("mult_e", 1.0f);
	AiParameterFlt("saturation", 1.0f);
	AiParameterBool("ray_reflected", 1);
	AiParameterBool("ray_refracted", 1);
	AiParameterBool("ray_diffuse", 1);
	AiParameterBool("ray_glossy", 1);
	AiParameterBool("ray_camera", 1);
	AiParameterFlt("candela_conversion", 100000.0f);
	AiParameterFlt("max_visibility_distance", 500.0f);
	AiParameterFlt("fog_transparency", 0.0f);
	AiParameterFlt("atmos_depth", 1.0f);
	AiParameterFlt("sun_inscatter", 1.0f);
	AiParameterFlt("sky_inscatter", 1.0f);
	AiParameterFlt("red_blue_shift", 0.0f);
	AiParameterFlt("tonemap", 0.0f);
	AiParameterFlt("inscatter_fade", 0.0f);
	AiParameterFlt("attenuation_fade", 0.0f);
	AiParameterBool("cie_overcast", 0);
	AiParameterBool("log_sun_intensity", 0);
	AiParameterFlt("camera_height", 1.0f);
}

#include "arnoldShaderFunctions.h"

node_initialize
{
	physicalSky *skyPtr = new physicalSky;
	AiNodeSetLocalData(node, skyPtr);
	AiMsgInfo("[aaPhysicalSky] Allocated sky data");
}

node_update
{
	physicalSky *skyPtr = (physicalSky *)AiNodeGetLocalData(node);
	AtVector v_sunDir;

	v_sunDir.x 		= AiNodeGetFlt(node, "sun_dir_x");
	v_sunDir.y 		= AiNodeGetFlt(node, "sun_dir_y");
	v_sunDir.z 		= AiNodeGetFlt(node, "sun_dir_z");
	bool _y_is_up 	= AiNodeGetBool(node, "y_is_up");
	Real _turbidity = AiNodeGetFlt(node, "turbidity");
	Real _multA 	= AiNodeGetFlt(node, "mult_a");
	Real _multB 	= AiNodeGetFlt(node, "mult_b");
	Real _multC 	= AiNodeGetFlt(node, "mult_c");
	Real _multD 	= AiNodeGetFlt(node, "mult_d");
	Real _multE 	= AiNodeGetFlt(node, "mult_e");
	Real _multiplier 	 = AiNodeGetFlt(node, "multiplier");
	Real _horizon_blur 	 = AiNodeGetFlt(node, "horizon_blur");
	Real _sun_opacity 	 = AiNodeGetFlt(node, "sun_opacity");
	Real _sun_intensity  = AiNodeGetFlt(node, "sun_intensity");
	Real _sun_disk_scale = AiNodeGetFlt(node, "sun_disk_scale");
	Real _saturation	 = AiNodeGetFlt(node, "saturation");
	Real _horizon_height = AiNodeGetFlt(node, "horizon_height");
	Real _candelaConversion = AiNodeGetFlt(node, "candela_conversion");
	Real _maxVisibilityDistance = AiNodeGetFlt(node, "max_visibility_distance");
	Real _fogTransparency = AiNodeGetFlt(node, "fog_transparency");
	Real _atmosDepth 	  = AiNodeGetFlt(node, "atmos_depth");
	Real _sunInscatter 	  = AiNodeGetFlt(node, "sun_inscatter");
	Real _skyInscatter 	  = AiNodeGetFlt(node, "sky_inscatter");
	Real _redBlueShift 	  = AiNodeGetFlt(node, "red_blue_shift");
	bool _cieOvercast 	  = AiNodeGetBool(node, "cie_overcast");
	AtRGBA _ground_color = AiNodeGetRGBA(node, "ground_color");
	
	
	if (v_sunDir.x == 0.0f && v_sunDir.y == 0.0f && v_sunDir.z == 0.0f)
		v_sunDir.y = 1.0;

	// scale sun vector to match horizon height GUI control
	v_sunDir = AiV3Normalize(v_sunDir);

	// z axis is up in spherical coords
	v_sunDir = cartesian_to_spherical_coords(v_sunDir);
	if (!_y_is_up)
	{
		Real y = v_sunDir.y;
		v_sunDir.y = v_sunDir.z;
		v_sunDir.z = y;
	}

	// Turbidity forced to start from 2.0
	// See "A Critical Review of the Preetham Skylight Model"
	// by Georg Zotti, Alexander Wilkie, Werner Purgathofer
	// http://wscg.zcu.cz/WSCG2007/Papers_2007/short/E59-full.pdf
	if (_turbidity < skyPtr->maxTurbidity)
		_turbidity = skyPtr->maxTurbidity;

	// attenuate sky horizon brightness to fix pink band near horizon when sun is high in sky
	//  does not work with Turbidity values lower than 2.0
	Real sunZ = v_sunDir.z;
	sunZ = std::clamp(sunZ, 0.0, 1.0);
	_multB = rescale(_multB, 0.0, 1.0, 0.45, 0.0, CLAMP_TO_RANGE);
	_multB = rescale(sunZ, 0.0, 1.0, 1.0, 0.55, CLAMP_TO_RANGE) + _multB;

	// shared variables access
	skyPtr->skyInput(_turbidity,
						v_sunDir,
						_multA,
						_multB,
						_multC,
						_multD,
						_multE,
						_ground_color,
						_multiplier,
						_horizon_height,
						_horizon_blur,
						_sun_opacity,
						_sun_intensity,
						_sun_disk_scale,
						_saturation,
						_candelaConversion,
						_maxVisibilityDistance,
						_fogTransparency,
						_atmosDepth,
						_sunInscatter,
						_skyInscatter,
						_redBlueShift,
						_cieOvercast);

		if(AiNodeGetBool(node, "log_sun_intensity"))
		{
			AtRGB sunColor = skyPtr->pSun.sunLightCol * skyPtr->candelaConversion * skyPtr->multiplier * skyPtr->sun_intensity;
			AiMsgWarning("[aaPhysicalSky] Sunlight color -- Red: %f, Green: %f, Blue: %f", sunColor.r, sunColor.g, sunColor.b);
		}

		far_clip = AiNodeGetFlt(AiUniverseGetCamera(AiNodeGetUniverse(node)), "far_clip");
		skyPtr->is_init = true;
}

shader_evaluate
{
	// evaluate ray visibility
	if (!isVisible(node, sg))
		return;

	physicalSky *skyPtr = (physicalSky *)AiNodeGetLocalData(node);

	// scale view vector to match horizon height GUI control
	AtVector direction = sg->Rd;
	direction.y += skyPtr->horizon_height;
	direction = AiV3Normalize(direction);

	// spherical coords setup
	Real thetav, phiv;
	direction = cartesian_to_spherical_coords(direction);
	if (!AiShaderEvalParamBool(p_y_is_up))
	{
		Real y = direction.y;
		direction.y = direction.z;
		direction.z = y;
	}
	skyPtr->computeAngles(phiv, thetav, direction);

	// Atmospheric Depth
	const Real s = 32000.0;
	// Height above ground, inverted for fog control
	const Real h0 = -1.0 * AiShaderEvalParamFlt(p_camera_height);

	// set up directions to fetch sky radiance
	AtVector skyDir = direction;
	skyDir.z = fabs(skyDir.z); // mirror sky

	// get sky radiance
	spectrum skyRadiance = skyPtr->getSkySpectralRadiance(skyDir);
	// get inscattered light in to the viewing direction
	spectrum inScatter = skyPtr->ap_InscatteredRadiance(h0, s, direction);
	// convert spectrums to RGB
	AtRGB inScatterRGB = inScatter.toRGB();
	inScatter *= (1.0 - skyPtr->fogTransparency);

	// Horizon attenuation
	// Don't apply attenuation spectrum if cieOvercast is On
	if (!AiShaderEvalParamBool(p_cie_overcast))
	{
		spectrum attenuation = skyPtr->ap_ExtinctionFactor(h0, s, direction);
		// change attenuation spectrum based on sun position
		attenuateSky(attenuation, skyPtr, (1.0 - AiShaderEvalParamFlt(p_attenuation_fade)));
		skyRadiance *= attenuation;
		skyRadiance += inScatter;
	}
	else
	{
		XYZ _XYZ = ACEScgtoXYZ(inScatterRGB);
		xyY _xyY = XYZtoxyY(_XYZ);
		_xyY.x = skyPtr->overcastWhitePoint_D65_x;
		_xyY.y = skyPtr->overcastWhitePoint_D65_x;
		xyYtoRGB(_xyY, inScatterRGB);
		skyRadiance += inScatter;
	}

	// convert spectrums to RGB
	AtRGB skyRadianceRGB = skyRadiance.toRGB();

	// fade in inScatter
	fadeInScatter(AiShaderEvalParamFlt(p_inscatter_fade), skyRadianceRGB, inScatterRGB);

	// apply transparency to inScatter
	AtRGB inScatterTRGB = inScatterRGB * (1.0 - skyPtr->fogTransparency);

	// add sun and ground
	skyPtr->addGround(skyRadianceRGB, inScatterTRGB, direction);
	skyPtr->addSun(skyRadianceRGB, inScatterTRGB, direction);

	// do the candela conversions, and also check for user-defined multiplier
	inScatterRGB *= skyPtr->candelaConversion * skyPtr->multiplier;
	skyRadianceRGB *= skyPtr->candelaConversion * skyPtr->multiplier;

	skyRadianceRGB = getVolumeEffects(skyPtr, skyRadianceRGB, inScatterRGB, sg);

	// fade sky based on sun's position below/near horizon
	// allow for about 6 degrees below horizon of light
	dusk(skyRadianceRGB, skyPtr);

	// apply saturation control
	skyRadianceRGB = applySaturation(skyRadianceRGB, skyPtr->saturation);

	// tonemap
	skyRadianceRGB = tonemap(skyRadianceRGB, AiShaderEvalParamFlt(p_tonemap));

	sg->out.RGB() = skyRadianceRGB;
	
}

node_finish
{
	delete (physicalSky *)AiNodeGetLocalData(node);
	AiMsgInfo("[aaPhysicalSky] Deallocated sky data");
}

node_loader
{
	if (i > 0)
		return false;

	node->methods = aaPhysicalSkyMethods;
	node->output_type = AI_TYPE_RGB;
	node->name = "aaPhysicalSky";
	node->node_type = AI_NODE_SHADER;
	strcpy(node->version, AI_VERSION);
	return true;
}
