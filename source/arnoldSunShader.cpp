//#include "ai.h"
//#include <string>
//#include "math.h"
//#include "physicalSunClass.h"
//
//AI_SHADER_NODE_EXPORT_METHODS(aaPhysicalSunMethods);
//
//enum aaPhysicalSunParams
//{
//	p_on,
//	p_turbidity,
//	p_sunDirX,
//	p_sunDirY,
//	p_sunDirZ,
//	p_multiplier,
//	p_exposure,
//};
//
//node_parameters
//{
//	AiParameterBOOL( "on"					, 1);
//	AiParameterFLT ( "turbidity"			, 5.0f);
//	AiParameterFLT ( "sunDirX"				, 0.0f);
//	AiParameterFLT ( "sunDirY"				, 0.0f);
//	AiParameterFLT ( "sunDirZ"				, 0.0f);
//	AiParameterFLT ( "multiplier"			, 1.0f);
//	AiParameterFLT ( "exposure"				, 1.0f);
//}
//
//node_initialize
//{
//}
//
//node_finish
//{
//}
//
//shader_evaluate
//{
//
//}
//
//node_loader
//{
//   if (i > 0)
//      return FALSE;
//
//   node->methods      = aaPhysicalSunMethods;
//   node->output_type  = AI_TYPE_RGB;
//   node->name         = "aaPhysicalSun";
//   node->node_type    = AI_NODE_SHADER;
//   strcpy(node->version, AI_VERSION);
//   return TRUE;
//}
