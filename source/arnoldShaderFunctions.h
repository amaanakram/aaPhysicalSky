#ifndef ARNOLDSHADERFUNCS
#define ARNOLDSHADERFUNCS

AtVector cartesian_to_spherical_coords(const AtVector &cart_rayDir)
{
	// assuming cart_rayDir is normalized
	AtVector v_spherical_coord_dir;
	v_spherical_coord_dir.y = cart_rayDir.z;
	v_spherical_coord_dir.z = cart_rayDir.y;
	v_spherical_coord_dir.x = cart_rayDir.x;
	return v_spherical_coord_dir;
}

bool isVisible(AtNode *&node, AtShaderGlobals *&sg)
{
	// ToDo: currently not working for sg->Ci

	bool result = 1;
	if (!AiShaderEvalParamBool(p_on))
		result = 0;

	// ray type checks for visibility
	if (sg->Rt == AI_RAY_CAMERA && !AiShaderEvalParamBool(p_ray_camera))
	{
		sg->out.RGBA() = AI_RGBA_ZERO;
		result = 0;
	}
	if (sg->Rt == AI_RAY_SPECULAR_REFLECT && !AiShaderEvalParamBool(p_ray_reflected))
		result = 0;
	if (sg->Rt == AI_RAY_SPECULAR_TRANSMIT && !AiShaderEvalParamBool(p_ray_refracted))
		result = 0;
	if (sg->Rt == AI_RAY_ALL_DIFFUSE && !AiShaderEvalParamBool(p_ray_diffuse))
		result = 0;
	if (sg->Rt == AI_RAY_ALL_SPECULAR && !AiShaderEvalParamBool(p_ray_glossy))
		result = 0;

	return result;
}

f_inline void attenuateSky(spectrum &attenuation, physicalSky *&skyPtr, Real user_t)
{
	Real t = RadsToDegs(skyPtr->thetaSun);
	t = rescale(t, 60.0, 90.0, 0.0, 1.0, CLAMP_TO_RANGE);
	clamp(user_t, 0.00001, 1.0); // not good
	t = minimum(t, user_t);
	clamp(t, 0.0, 1.0);
	attenuation.lerp(t, skyPtr->specOne);
}

f_inline void dusk(AtRGB &skyRadianceRGB, physicalSky *&skyPtr)
{
	Real t = RadsToDegs(skyPtr->thetaSun); // do all in radians eventually, and move to node_update
	t = rescale(t, 90.0, 96.0, 1.0, 0.0, CLAMP_TO_RANGE);
	skyRadianceRGB *= t;
}

AtRGB getVolumeEffects(physicalSky *&skyPtr, AtRGB &skyRadianceRGB, AtRGB &volumeRGB, AtShaderGlobals *&sg)
{
	AtRGB result;
	AtRGB unoccluded_color = AtRGB(1.0f, 0.0f, 0.2f); // sg->Ci;
	if (sg->sc == AI_CONTEXT_VOLUME)
	{
		// check for no-hit
		if (sg->Rl >= AI_BIG || sg->Rl > far_clip)
			result = skyRadianceRGB;
		// check for completely transparent fog
		else if (skyPtr->fogTransparency >= 1.0f)
			result = unoccluded_color;
		else
		{
			Real rayl = sg->Rl;
			// clamp raylength to be up to max visibility distance
			clamp(rayl, 0.0, skyPtr->maxVisibilityDistance);
			rayl = rescale(rayl, 0.0, skyPtr->maxVisibilityDistance, 0.0, 1.0, CLAMP_TO_RANGE);

			// fog out completely beyond maxVisibilityDistance
			volumeRGB = AiLerp(rayl, unoccluded_color, volumeRGB);

			// Transparency: when completey transparent (1.0), return unoccluded_color, if 0.0, return inScatter, lerp otherwise
			result = AiLerp(skyPtr->fogTransparency, volumeRGB, unoccluded_color);
		}
	}
	else
		result = skyRadianceRGB;

	return result;
}

inline Real tonemap_helper(Real x)
{
	// www.filmicgames.com/archives/75
	const float A = 0.15;
	const float B = 0.50;
	const float C = 0.10;
	const float D = 0.20;
	const float E = 0.02;
	const float F = 0.30;
	return ((x * (A * x + C * B) + D * E) / (x * (A * x + B) + D * F)) - E / F;
}

f_inline AtRGB tonemap(AtRGB &inColor, Real t)
{
	if (t == 0.0f)
		return inColor;
	else
	{
		XYZ _XYZ = ACEScgtoXYZ(inColor);
		xyY _xyY = XYZtoxyY(_XYZ);

		float ExposureBias = 2.0f;
		float curr = tonemap_helper(ExposureBias * _xyY.Y);

		const float W = 11.2;
		float whiteScale = 1.0f / tonemap_helper(W);
		_xyY.Y = curr * whiteScale;

		AtRGB result;
		xyYtoRGB(_xyY, result);

		clamp(t, 0.0, 1.0);
		result = AiLerp(t, inColor, result);
		return result;
	}
}

f_inline AtRGB applySaturation(AtRGB &in, Real saturation)
{
	if (saturation == 1.0)
		return in;

	AtRGB result;
	Real r, g, b;
	Real h, s, v;
	RGBtoHSV(in.r, in.g, in.b, &h, &s, &v);
	s *= saturation;
	HSVtoRGB(&r, &g, &b, h, s, v);
	result.r = r;
	result.g = g;
	result.b = b;
	return result;
}

f_inline void fadeInScatter(Real fade, AtRGB &skyRadianceRGB, AtRGB &inScatterRGB)
{
	if (fade > 0.0f)
	{
		clamp(fade, 0.0, 1.0);
		skyRadianceRGB = AiLerp(fade, skyRadianceRGB, inScatterRGB);
	}
	if (fade == 1.0f)
		skyRadianceRGB = inScatterRGB;
}
#endif // ARNOLDSHADERFUNCS_H