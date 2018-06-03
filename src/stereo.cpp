#include "stereo.h"
#include "tables.h"




//
//impl

#define C_INV_SQRT_2 0.70710678118654752440

namespace OpenMP3
{

	void StereoIntensityLong(UInt sfreq, const FrameData & data, UInt gr, UInt sfb, Float32 is[2][576]);

	void StereoIntensityShort(UInt sfreq, const FrameData & data, UInt gr, UInt sfb, Float32 is[2][576]);

	extern const float kIsRatios[6];

}




//
//

void OpenMP3::Stereo(UInt sfreq, UInt8 joint_stereo_mode, const FrameData & data, unsigned gr, Float32 is[2][576])
{
	//Do nothing if joint stereo is not enabled
	//if ((header.mode != 1) || (header.mode_extension == 0)) return;

	if (joint_stereo_mode & 0x2) //Mid/Side stereo processing
	{
		/* Determine how many frequency lines to transform */
		UInt max_pos = data.count1[gr][!(data.count1[gr][0] > data.count1[gr][1])];
	
		/* Do the actual processing */
		for (UInt i = 0; i < max_pos; i++)
		{
			Float32 left = Float32((is[0][i] + is[1][i]) * C_INV_SQRT_2);
			
			Float32 right = Float32((is[0][i] - is[1][i]) * C_INV_SQRT_2);
			
			is[0][i] = left;
			
			is[1][i] = right;
		}
	}

	if (joint_stereo_mode & 0x1) //intensity stereo processing
	{
		//unsigned sfreq = header.sampling_frequency;

		//First band that is intensity stereo encoded is first band scale factor band on or above count1 frequency line. 
		//N.B.: Intensity stereo coding is only done for higher subbands, but logic is here for lower subbands
		if (data.window_switching[gr][0] && (data.block_type[gr][0] == 2))
		{ 
			//Short blocks
			
			//Check if the first two subbands *(=2*18 samples = 8 long or 3 short sfb's) uses long blocks
			if (data.mixed_block[gr][0])
			{
				UInt sfb;

				//First process 8 sfb's at start, where scale factor band above count1 for the right channel
				for (sfb = 0; sfb < 8; sfb++) if (kScaleFactorBandIndices[sfreq].l[sfb] >= data.count1[gr][1]) StereoIntensityLong(sfreq, data, gr, sfb, is);
				 
				//the remaining bands which uses short blocks
				for (sfb = 3; sfb < 12; sfb++) if (kScaleFactorBandIndices[sfreq].s[sfb] * 3 >= data.count1[gr][1]) StereoIntensityShort(sfreq, data, gr, sfb, is);
			}
			else
			{ 
				//Only short blocks, where scale factor band above count1 for the right channel
				for (UInt sfb = 0; sfb < 12; sfb++) if (kScaleFactorBandIndices[sfreq].s[sfb] * 3 >= data.count1[gr][1]) StereoIntensityShort(sfreq, data, gr, sfb, is);
			}
		}
		else 
		{
			//Only long blocks, where scale factor band above count1 for the right channel
			for (UInt sfb = 0; sfb < 21; sfb++) if (kScaleFactorBandIndices[sfreq].l[sfb] >= data.count1[gr][1]) StereoIntensityLong(sfreq, data, gr, sfb, is);
		}
	}
}

//intensity stereo processing for entire subband with long blocks
void OpenMP3::StereoIntensityLong(UInt sfreq, const FrameData & data, unsigned gr, unsigned sfb, Float32 is[2][576])
{
	unsigned i, sfb_start, sfb_stop, is_pos;
	float is_ratio_l, is_ratio_r;

	/* Check that((is_pos[sfb]=scalefac) != 7) => no intensity stereo */
	if ((is_pos = data.scalefac_l[gr][0][sfb]) != 7) 
	{
		//sfreq = g_frame_header.sampling_frequency; /* Setup sampling freq index */

		sfb_start = kScaleFactorBandIndices[sfreq].l[sfb];
		
		sfb_stop = kScaleFactorBandIndices[sfreq].l[sfb + 1];

		if (is_pos == 6) 
		{ /* tan((6*PI)/12 = PI/2) needs special treatment! */
			is_ratio_l = 1.0f;
			is_ratio_r = 0.0f;
		}
		else 
		{
			is_ratio_l = kIsRatios[is_pos] / (1.0f + kIsRatios[is_pos]);
			
			is_ratio_r = 1.0f / (1.0f + kIsRatios[is_pos]);
		}

		/* Now decode all samples in this scale factor band */
		for (i = sfb_start; i < sfb_stop; i++) 
		{
			float left = is_ratio_l * data.is[gr][0][i];

			float right = is_ratio_r * data.is[gr][0][i];
			
			is[0][i] = left;

			is[1][i] = right;
		}
	}
}

//intensity stereo processing for an entire subband that uses short blocks
void OpenMP3::StereoIntensityShort(UInt sfreq, const FrameData & data, unsigned gr, unsigned sfb, Float32 is[2][576])
{
	//TODO this code does seem to do anything, except for rounding the is buffer!!

	UInt is_pos, is_ratio_l, is_ratio_r;

	UInt win_len = kScaleFactorBandIndices[sfreq].s[sfb + 1] - kScaleFactorBandIndices[sfreq].s[sfb];
	
	/* The three windows within the band has different scalefactors */
	for (unsigned win = 0; win < 3; win++)
	{
		/* Check that((is_pos[sfb]=scalefac) != 7) => no intensity stereo */
		if ((is_pos = data.scalefac_s[gr][0][sfb][win]) != 7) 
		{
			UInt sfb_start = kScaleFactorBandIndices[sfreq].s[sfb] * 3 + win_len * win;

			UInt sfb_stop = sfb_start + win_len;

			if (is_pos == 6) //tan((6*PI)/12 = PI/2) needs special treatment 
			{ 
				is_ratio_l = 1;

				is_ratio_r = 0;
			}
			else 
			{
				is_ratio_l = UInt(kIsRatios[is_pos] / (1.0f + kIsRatios[is_pos]));

				is_ratio_r = UInt(1.0f / (1.0f + kIsRatios[is_pos]));
			}
			
			/* Now decode all samples in this scale factor band */
			for (UInt i = sfb_start; i < sfb_stop; i++)
			{
				//Float32 left = is_ratio_l = data.is[gr][0][i];

				//Float32 right = is_ratio_r = data.is[gr][0][i];

				Float32 left = Float32(is_ratio_l = UInt(data.is[gr][0][i]));

				Float32 right = Float32(is_ratio_r = UInt(data.is[gr][0][i]));

				is[0][i] = left;

				is[1][i] = right;
			}
		} 
	} 
}




//
//

const float OpenMP3::kIsRatios[6] = { 0.000000f,0.267949f,0.577350f,1.000000f,1.732051f,3.732051f };
