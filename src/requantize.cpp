#include "requantize.h"

#include "types.h"
#include "tables.h"




//
//internal

namespace OpenMP3
{

	Float32 Requantize_Pow_43(Float32 is_pos);

	void RequantizeLong(const FrameData & data, UInt gr, UInt ch, UInt is_pos, UInt sfb, Float32 is[576]);

	void RequantizeShort(const FrameData & data, UInt gr, UInt ch, UInt is_pos, UInt sfb, UInt win, Float32 is[576]);


	extern const Float32 kPreTab[21];

}




//
//impl

void OpenMP3::Requantize(UInt sfreq, const FrameData & data, UInt gr, UInt ch, Float32 is[576])
{
	UInt sfb /* scalefac band index */, next_sfb /* frequency of next sfb */, i, j;

	/* Determine type of block to process */
	if (data.window_switching[gr][ch] && (data.block_type[gr][ch] == 2))	//Short blocks
	{ 
		/* Check if the first two subbands (=2*18 samples = 8 long or 3 short sfb's) uses long blocks */
		if (data.mixed_block[gr][ch])
		{ 
			/* 2 longbl. sb  first */
			/* First process the 2 long block subbands at the start */

			sfb = 0;

			next_sfb = kScaleFactorBandIndices[sfreq].l[sfb + 1];

			for (i = 0; i < 36; i++) 
			{
				if (i == next_sfb) 
				{
					sfb++;

					next_sfb = kScaleFactorBandIndices[sfreq].l[sfb + 1];
				}

				RequantizeLong(data, gr, ch, i, sfb, is);
			}
		
			/* And next the remaining,non-zero,bands which uses short blocks */
			sfb = 3;
			
			next_sfb = kScaleFactorBandIndices[sfreq].s[sfb + 1] * 3;
			
			UInt win_len = kScaleFactorBandIndices[sfreq].s[sfb + 1] - kScaleFactorBandIndices[sfreq].s[sfb];

			for (i = 36; i < data.count1[gr][ch]; /* i++ done below! */)
			{
				/* Check if we're into the next scalefac band */
				if (i == next_sfb)
				{
					ASSERT(sfb < 14);

					sfb++;

					next_sfb = kScaleFactorBandIndices[sfreq].s[sfb + 1] * 3;

					win_len = kScaleFactorBandIndices[sfreq].s[sfb + 1] - kScaleFactorBandIndices[sfreq].s[sfb];
				} 

				for (UInt win = 0; win < 3; win++) for (j = 0; j < win_len; j++) RequantizeShort(data, gr, ch, i++, sfb, win, is);
			}
		}
		else //Only short blocks
		{ 
			sfb = 0;

			next_sfb = kScaleFactorBandIndices[sfreq].s[sfb + 1] * 3;

			UInt win_len = kScaleFactorBandIndices[sfreq].s[sfb + 1] - kScaleFactorBandIndices[sfreq].s[sfb];

			for (i = 0; i < data.count1[gr][ch]; /* i++ done below! */)
			{
				/* Check if we're into the next scalefac band */
				if (i == next_sfb) 
				{
					ASSERT(sfb < 14);

					sfb++;
					
					next_sfb = kScaleFactorBandIndices[sfreq].s[sfb + 1] * 3;
					
					win_len = kScaleFactorBandIndices[sfreq].s[sfb + 1] - kScaleFactorBandIndices[sfreq].s[sfb];
				}

				for (UInt win = 0; win < 3; win++) for (j = 0; j < win_len; j++) RequantizeShort(data, gr, ch, i++, sfb, win, is);
			}
		}
	}
	else //Only long blocks
	{ 
		sfb = 0;

		next_sfb = kScaleFactorBandIndices[sfreq].l[sfb + 1];

		UInt n = data.count1[gr][ch];

		for (i = 0; i < n; i++) 
		{
			if (i == next_sfb) 
			{
				ASSERT(sfb < 23);

				sfb++;

				next_sfb = kScaleFactorBandIndices[sfreq].l[sfb + 1];
			}

			RequantizeLong(data, gr, ch, i, sfb, is);
		}
	}
}


void OpenMP3::Reorder(UInt sfreq, const FrameData & data, UInt gr, UInt ch, Float32 is[576])
{
	OpenMP3::Float32 re[576];	//TODO use working buffer

	if (data.window_switching[gr][ch] && (data.block_type[gr][ch] == 2))	//Only reorder short blocks
	{
		/* Check if the first two subbands (=2*18 samples = 8 long or 3 short sfb's) uses long blocks */

		auto s = kScaleFactorBandIndices[sfreq].s;

		UInt sfb = data.mixed_block[gr][ch] ? 3 : 0; /* 2 longbl. sb  first */

		UInt count1 = data.count1[gr][ch];

		UInt next_sfb = s[sfb + 1] * 3;

		UInt win_len = s[sfb + 1] - s[sfb];

		for (UInt i = ((sfb == 0) ? 0 : 36); i < 576; /* i++ done below! */)
		{
			/* Check if we're into the next scalefac band */
			if (i == next_sfb)
			{
				/* Copy reordered data back to the original vector */
				for (UInt j = 0; j < 3 * win_len; j++) is[3 * s[sfb] + j] = re[j];

				/* Check if this band is above the rzero region,if so we're done */
				if (i >= count1) return;

				sfb++;

				next_sfb = s[sfb + 1] * 3;

				win_len = s[sfb + 1] - s[sfb];
			}

			//Do the actual reordering

			for (UInt win = 0; win < 3; win++) for (UInt j = 0; j < win_len; j++) re[j * 3 + win] = is[i++];
		}

		/* Copy reordered data of last band back to original vector */
		for (UInt j = 0; j < 3 * win_len; j++) is[3 * s[12] + j] = re[j];
	}
}

inline float OpenMP3::Requantize_Pow_43(Float32 f_is_pos)
{
	int is_pos = int(f_is_pos);

	return powf((float)is_pos, 4.0f / 3.0f);
}

//requantize sample in subband that uses long blocks
void OpenMP3::RequantizeLong(const FrameData & data, UInt gr, UInt ch, UInt is_pos, UInt sfb, Float32 is[576])
{
	ASSERT(is_pos < 576);

	float sf_mult = data.scalefac_scale[gr][ch] ? 1.0f : 0.5f;
	
	float pf_x_pt = data.preflag[gr][ch] * kPreTab[sfb];
	
	Float64 tmp1 = pow(2.0, -(sf_mult *(data.scalefac_l[gr][ch][sfb] + pf_x_pt)));
	
	Float64 tmp2 = pow(2.0, 0.25 * (data.global_gain[gr][ch] - 210.0));

	Float64 tmp3;

	if (is[is_pos] < 0.0)
	{
		tmp3 = -Requantize_Pow_43(-is[is_pos]);
	}
	else
	{
		tmp3 = Requantize_Pow_43(is[is_pos]);
	}

	is[is_pos] = Float32(tmp1 * tmp2 * tmp3);
}

//requantize sample in subband that uses short blocks
void OpenMP3::RequantizeShort(const FrameData & data, UInt gr, UInt ch, UInt is_pos, UInt sfb, UInt win, Float32 is[576])
{
	ASSERT(is_pos < 576);

	Float32 sf_mult = data.scalefac_scale[gr][ch] ? 1.0f : 0.5f;

	Float64 tmp1 = pow(2.0, -(sf_mult * data.scalefac_s[gr][ch][sfb][win]));

	Float64 tmp2 = pow(2.0, 0.25 * (data.global_gain[gr][ch] - 210.0f - 8.0f * data.subblock_gain[gr][ch][win]));

	Float64 tmp3 = (is[is_pos] < 0.0) ? -Requantize_Pow_43(-is[is_pos]) : Requantize_Pow_43(is[is_pos]);

	is[is_pos] = Float32(tmp1 * tmp2 * tmp3);
}




//
//

const float OpenMP3::kPreTab[21] = { 0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,2,3,3,3,2 };
