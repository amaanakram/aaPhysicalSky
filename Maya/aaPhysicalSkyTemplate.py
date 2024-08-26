import pymel.core as pm
import mtoa.utils as utils
import mtoa.ui.ae.utils as aeUtils
from mtoa.ui.ae.shaderTemplate import ShaderAETemplate

class AEaaPhysicalSkyTemplate(ShaderAETemplate):
	
	def setup(self):
		self.addSwatch()

		self.beginScrollLayout()

		self.addCustom('message', 'AEshaderTypeNew','AEshaderTypeReplace')

		self.beginLayout("Global Parameters", collapse=False)
		self.addControl("on", label="On")
		self.addControl("y_is_up", label="Y is UP")
		self.addControl("multiplier", label="Multiplier")
		self.addControl("saturation", label="Saturation")
		self.addControl("redBlueShift", label="Red-Blue Shift")
		self.addControl("multB", label="Horizon Luma")
		self.addControl("tonemap", label="Tone Map")
		self.endLayout()
		
		self.beginLayout("Atmospheric and Volume Parameters", collapse=True)
		self.addControl("turbidity", label="Haze")
		self.addControl("horizon_height", label="Horizon Height")
		self.addControl("horizon_blur", label="Horizon Blur")
		self.addControl("ground_color", label="Ground Colour")
		self.addControl("maxVisibilityDistance", label="Max Visibility Distance")
		self.addControl("fogTransparency", label="Transparency")
		self.addControl("cameraHeight", label="Camera Height")
		self.addControl("sunInscatter", label="Sun Scattering")
		self.addControl("skyInscatter", label="Sky Scattering")
		self.endLayout()
		
		self.beginLayout("Sun Parameters", collapse=True)
		self.addControl("sunDirX", label="Direction X")
		self.addControl("sunDirY", label="Direction Y")
		self.addControl("sunDirZ", label="Direction Z")
		self.addControl("sun_intensity", label="Intensity")
		self.addControl("sun_opacity", label="Opacity")
		self.addControl("sun_disk_scale", label="Size")
		self.addControl("multD", label="Glow")
		self.endLayout()
		
		self.beginLayout("Visibility", collapse=True)
		self.addControl("rayCamera", label="Camera Rays")
		self.addControl("rayReflected", label="Reflected Rays")
		self.addControl("rayRefracted", label="Refracted Rays")
		self.addControl("rayDiffuse", label="Diffuse Rays")
		self.addControl("rayGlossy", label="Glossy Rays")
		self.endLayout()
		
		self.beginLayout("Misc", collapse=True)
		self.addControl("inScatterFade", label="Fade in In-Scatter")
		self.addControl("attenuationFade", label="Fade out Attenuation")
		self.addControl("candelaConversion", label="Cd/m^2 Conversion")
		self.addControl("cieOvercast", label="CIE Overcast Sky Model")
		self.addControl("logSunIntensity", label="Log Sun Intensity in Console")
		self.endLayout()

		pm.mel.AEdependNodeTemplate(self.nodeName)
