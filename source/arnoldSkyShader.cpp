#include <string.h>
#include "ai.h"
#include "physicalSkyClass.h"


static float far_clip; // camera far clipping distance

AI_SHADER_NODE_EXPORT_METHODS(aaPhysicalSkyMethods);

enum aaPhysicalSkyParams
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
    AiParameterBOOL( "on"                   , 1);
    AiParameterBOOL( "y_is_up"              , 1);
	AiParameterFLT ( "turbidity"			, 5.0f);
	AiParameterFLT ( "sun_dir_x"			, 0.0f);
	AiParameterFLT ( "sun_dir_y"			, 0.0f);
	AiParameterFLT ( "sun_dir_z"			, 0.0f);
	AiParameterFLT ( "multiplier"			, 1.0f);
	AiParameterRGBA( "ground_color"			, 0.2f, 0.2f, 0.2f, 1.0f);
	AiParameterFLT ( "horizon_height"		, 0.0f);
	AiParameterFLT ( "horizon_blur"			, 0.1f);
	AiParameterFLT ( "sun_opacity"			, 1.0f);
	AiParameterFLT ( "sun_disk_scale"		, 1.0f);
	AiParameterFLT ( "sun_intensity"		, 10.0f);
	AiParameterFLT ( "mult_a"				, 1.0f);
	AiParameterFLT ( "mult_b"				, 1.0f);
	AiParameterFLT ( "mult_c"				, 1.0f);
	AiParameterFLT ( "mult_d"				, 1.0f);
	AiParameterFLT ( "mult_e"				, 1.0f);
	AiParameterFLT ( "saturation"			, 1.0f);
	AiParameterBOOL( "ray_reflected"		, 1);
	AiParameterBOOL( "ray_refracted"		, 1);
	AiParameterBOOL( "ray_diffuse"			, 1);
	AiParameterBOOL( "ray_glossy"			, 1);
	AiParameterBOOL( "ray_camera"			, 1);
	AiParameterFLT ( "candela_conversion"	, 100000.0f);
	AiParameterFLT ( "max_visibility_distance", 500.0f);
	AiParameterFLT ( "fog_transparency"		, 0.0f);
	AiParameterFLT ( "atmos_depth"			, 1.0f);
	AiParameterFLT ( "sun_inscatter"		, 1.0f);
	AiParameterFLT ( "sky_inscatter"		, 1.0f);
	AiParameterFLT ( "red_blue_shift"		, 0.0f);
	AiParameterFLT ( "tonemap"				, 0.0f);
	AiParameterFLT ( "inscatter_fade"		, 0.0f);
	AiParameterFLT ( "attenuation_fade"		, 0.0f);
	AiParameterBOOL( "cie_overcast"			, 0);
	AiParameterBOOL( "log_sun_intensity"	, 0);
	AiParameterFLT ( "camera_height"		, 1.0f);
}

#include "arnoldShaderFunctions.h"

node_initialize
{
	physicalSky *skyPtr = new physicalSky;
	AiNodeSetLocalData(node,skyPtr);
	AiMsgInfo("[aaPhysicalSky] Allocated sky data");
}

node_update
{
	physicalSky *skyPtr = (physicalSky *)AiNodeGetLocalData(node);

    bool    _y_is_up  = params[p_y_is_up].BOOL;
	Real	_turbidity =  params[p_turbidity].FLT;
	Real	_multA = params[p_mult_a].FLT;
    Real	_multB = params[p_mult_b].FLT;
    Real	_multC = params[p_mult_c].FLT;
    Real	_multD = params[p_mult_d].FLT;
    Real	_multE = params[p_mult_e].FLT;
	AtRGBA	_ground_color = params[p_ground_color].RGBA;
    Real	_multiplier = params[p_multiplier].FLT;
    Real	_horizon_blur = params[p_horizon_blur].FLT;
    Real	_sun_opacity = params[p_sun_opacity].FLT;
    Real	_sun_intensity = params[p_sun_intensity].FLT;
    Real	_sun_disk_scale = params[p_sun_disk_scale].FLT;
    Real	_saturation = params[p_saturation].FLT;
	Real	_horizon_height = params[p_horizon_height].FLT;
	Real	_candelaConversion = params[p_candela_conversion].FLT;
	Real	_maxVisibilityDistance = params[p_max_visibility_distance].FLT;
	Real	_fogTransparency = params[p_fog_transparency].FLT;
	Real	_atmosDepth = params[p_atmos_depth].FLT;
	Real	_sunInscatter = params[p_sun_inscatter].FLT;
	Real	_skyInscatter = params[p_sky_inscatter].FLT;
	Real	_redBlueShift = params[p_red_blue_shift].FLT;
	bool	_cieOvercast  = params[p_cie_overcast].BOOL;

	AtVector v_sunDir;
	v_sunDir.x	= params[p_sun_dir_x].FLT;
	v_sunDir.y	= params[p_sun_dir_y].FLT;
	v_sunDir.z	= params[p_sun_dir_z].FLT;
	if(v_sunDir.x == 0.0f && v_sunDir.y == 0.0f && v_sunDir.z == 0.0f)
		v_sunDir.y = 1.0;

	// scale sun vector to match horizon height GUI control
	AiV3Normalize(v_sunDir,v_sunDir);

	//z axis is up in spherical coords
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
	if(_turbidity < skyPtr->maxTurbidity)
		_turbidity = skyPtr->maxTurbidity;

	 //attenuate sky horizon brightness to fix pink band near horizon when sun is high in sky
	// does not work with Turbidity values lower than 2.0
	Real sunZ = v_sunDir.z;
	clamp(sunZ, 0.0, 1.0);
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

	if(params[p_log_sun_intensity].BOOL)
	{
		AtColor sunColor = skyPtr->pSun.sunLightCol * skyPtr->candelaConversion * skyPtr->multiplier * skyPtr->sun_intensity;
		AiMsgWarning("[aaPhysicalSky] Sunlight color -- Red: %f, Green: %f, Blue: %f", sunColor.r, sunColor.g, sunColor.b);
	}

	far_clip = AiNodeGetFlt(AiUniverseGetCamera(), "far_clip");
}

shader_evaluate
{
	// evaluate ray visibility
	if(!isVisible(node, sg))
		return;

	physicalSky *skyPtr = (physicalSky *)AiNodeGetLocalData(node);

	// scale view vector to match horizon height GUI control
	AtVector direction = sg->Rd;
	direction.y += skyPtr->horizon_height;
	AiV3Normalize(direction, direction);

	// spherical coords setup
	Real thetav, phiv;
	direction = cartesian_to_spherical_coords(direction);
	if (!AiShaderEvalParamBool(p_y_is_up))
	{
	  Real y = direction.y;
	  direction.y = direction.z;
	  direction.z = y;
	}
	skyPtr->computeAngles(phiv,thetav,direction);

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
	AtColor inScatterRGB = inScatter.toRGB();
	inScatter *= (1.0 - skyPtr->fogTransparency);
	
	// Horizon attenuation
	// Don't apply attenuation spectrum if cieOvercast is On
	if(!AiShaderEvalParamBool(p_cie_overcast))
	{
		spectrum attenuation = skyPtr->ap_ExtinctionFactor(h0, s, direction);
		// change attenuation spectrum based on sun position
		attenuateSky(attenuation, skyPtr, (1.0 - AiShaderEvalParamFlt(p_attenuation_fade)));
		skyRadiance	*= attenuation;
		skyRadiance += inScatter;
	}
	else
	{
		XYZ _XYZ = RGBtoXYZ(inScatterRGB);
		xyY _xyY = XYZtoxyY(_XYZ);
		_xyY.x = skyPtr->overcastWhitePoint_D65_x;
		_xyY.y = skyPtr->overcastWhitePoint_D65_x;
		xyYtoRGB(_xyY, inScatterRGB);
		skyRadiance	+= inScatter;
	}

	// convert spectrums to RGB
	AtColor skyRadianceRGB = skyRadiance.toRGB();

	// fade in inScatter
	fadeInScatter(AiShaderEvalParamFlt(p_inscatter_fade), skyRadianceRGB, inScatterRGB);

	// apply transparency to inScatter
	AtColor inScatterTRGB	= inScatterRGB * (1.0 - skyPtr->fogTransparency);
	
	// add sun and ground
	skyPtr->addGround(skyRadianceRGB, inScatterTRGB, direction);
	skyPtr->addSun(skyRadianceRGB, inScatterTRGB, direction);

	// do the candela conversions, and also check for user-defined multiplier
	inScatterRGB   *= skyPtr->candelaConversion * skyPtr->multiplier;
	skyRadianceRGB *= skyPtr->candelaConversion * skyPtr->multiplier;

	skyRadianceRGB = getVolumeEffects(skyPtr, skyRadianceRGB, inScatterRGB, sg);

	// fade sky based on sun's position below/near horizon
	// allow for about 6 degrees below horizon of light
	dusk(skyRadianceRGB, skyPtr);

	//apply saturation control
	skyRadianceRGB = applySaturation(skyRadianceRGB, skyPtr->saturation);

	// tonemap
	skyRadianceRGB = tonemap(skyRadianceRGB, AiShaderEvalParamFlt(p_tonemap));

	sg->out.RGB = skyRadianceRGB;

}

node_finish
{
	delete ( physicalSky*)AiNodeGetLocalData(node);
	AiMsgInfo("[aaPhysicalSky] Deallocated sky data");
}

node_loader
{
   if (i > 0)
      return FALSE;

   node->methods      = aaPhysicalSkyMethods;
   node->output_type  = AI_TYPE_RGB;
   node->name         = "aaPhysicalSky";
   node->node_type    = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
}
