#ifndef SPECTRUM_GLOBALS
#define SPECTRUM_GLOBALS

static spectrum netaTable[1801];

static Real netaLambdaTable[wavelengths][1801] =
{
	#include "netaTable.h"
};

static void initNetaTable()
{
    Real data[wavelengths];
	#pragma omp parallel for private(data)
	for(int i = 0; i < 1801; i++) 
	{
		for(int j = 0; j < wavelengths; j++) 
			data[j] = netaLambdaTable[j][i];
		netaTable[i] = spectrum(data); 
	}
}

#endif // SPECTRUM_GLOBALS_H