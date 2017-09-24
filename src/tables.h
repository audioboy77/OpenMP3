#pragma once




//
//declarations

namespace OpenMP3
{

	struct ScaleFactorBandIndices		//Scale factor band indices,for long and short windows
	{
		UInt l[23];
		UInt s[14];
	};

	
	extern const UInt kBitRates[15];

	extern const UInt kSampleRates[3];


	extern const ScaleFactorBandIndices kScaleFactorBandIndices[3];

	extern const UInt kScaleFactorSizes[16][2];


	extern const UInt8 kInfo[4];

}