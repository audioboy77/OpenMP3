#pragma once

#include "types.h"




//
//decl

namespace OpenMP3
{

	void Antialias(const FrameData & data, UInt gr, UInt ch, Float32 is[576]);

	void HybridSynthesis(const FrameData & data, UInt gr, UInt ch, Float32 store[2][32][18], Float32 is[576]);

	void FrequencyInversion(Float32 is[576]);

	void SubbandSynthesis(const FrameData & data, UInt ch, const Float32 is[576], Float32 v_vec[2][1024], Float32 output[576]);

}