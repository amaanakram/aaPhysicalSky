#ifndef CONSTANTS
#define CONSTANTS

/////////////////// TYPEDEFS ////////////////////////
#ifndef Real
#define Real double // switching to double-precision because of precision issues with floats and paper data
#endif
/* ----------------------------------------------- */

///////////////// HELPER DEFINES /////////////////
#define EXP(x) t2exp(x) //faster exponental function than the standard exp(x), error rate seems acceptable

#define CLAMP_TO_RANGE 1

#ifndef FLT_MIN
#define FLT_MIN 1.175494351e-38F /* min positive value */
#endif
/* ----------------------------------------------- */

/////////// STANDARD CONSTANTS /////////////////
static const Real aa_EPSILON = FLT_MIN;
static const Real aa_PI		 = 3.14159265358979323846;
static const Real aa_TWOPI	 = 6.28318530717958647692;
static const Real aa_PIBYTWO = 1.57079632679489661923;
static const Real aa_180BYPI = 57.295779513082320876798154814105;
static const Real aa_PIBY180 = 0.01745329251994329576922222222222;
static const Real aa_INV_PIBYTWO = 0.63661977236758134307607071493546;
/* ----------------------------------------------------------------- */

/////////////////////WAVELENGTH HELPER DEFINES///////////////////////////////

static const int nTheta = 20;				// Number of bins for theta
static const int nPhi	= 20;				// Number of bins for phi
static const int wavelengths	 = 38;		// number to wavelengths to process, 380nm to 750nm
static const int wavelengthDelta = 10;		// spacing of 10nm between wavelengths
/* ------------------------------------------------------------------------------ */

///////////// SSE/AVX intrinsics defines //////////////////////
#define BOUNDARY 16 // alignment boundary

#ifdef _MSC_VER
#define ALIGN(x) __declspec(align(x))
#else // gcc
#define ALIGN(x) __attribute__((aligned(x)))
#endif

#ifdef USE_SSE2
#include <emmintrin.h>
typedef __m128d v2d;
#endif

#ifdef USE_AVX
// Header for AVX intrinsics (includes SSE intrinsics.)
#include <immintrin.h>
typedef __m128d v2d;
#endif
/* ----------------------------------------------- */

//////////////// FORCE INLINE //////////////////////
#ifdef _MSC_VER
#define f_inline __forceinline
#else // gcc
#define f_inline inline
#endif
/* ----------------------------------------------- */

#endif // CONSTANTS_H