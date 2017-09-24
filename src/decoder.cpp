#include "../include/openmp3.h"

#include "types.h"
#include "tables.h"
#include "requantize.h"
#include "stereo.h"
#include "synthesis.h"
#include "huffman.h"





//nsa debug

#define C_PI 3.14159265358979323846




//
//private

struct OpenMP3::Decoder::Private
{
	static bool ReadSideInfo(Decoder & self, const Frame & header);

	static bool ReadMain(Decoder & self, const Frame & header);

	static void DecodeFrame(Decoder & self, const Frame & header, Float32 is[2][2][576], Float32 out576[2][1152]);

	static void ReadBytes(Decoder & self, UInt no_of_bytes, UInt data_vec[]);

	static const UInt8 * ReadBytes(Decoder & self, UInt nyytes);
};




//
//

OpenMP3::Library::Library()
{
	//for (UInt i = 0; i < 64; i++) for (UInt j = 0; j < 32; j++) m_sbs_n_win[i][j] = Float32(cos(Float64((16 + i) * (2 * j + 1)) * (C_PI / 64.0)));
}




//
//impl

OpenMP3::Decoder::Decoder(const Library & library)
	: library(library)
{
	CT_ASSERT(sizeof(m_br) >= sizeof(Reservoir));

	CT_ASSERT(sizeof(m_framedata) >= sizeof(FrameData));

	Reset();
}

void OpenMP3::Decoder::Reset()
{
	MemClear(m_br, sizeof(Reservoir));

	MemClear(m_framedata, sizeof(FrameData));

	MemClear(m_hs_store, sizeof(m_hs_store));

	MemClear(m_sbs_v_vec, sizeof(m_sbs_v_vec));

	//hsynth_init = synth_init = 1;
}

OpenMP3::UInt OpenMP3::Decoder::ProcessFrame(const Frame & frame, Float32 output[2][1152])
{
	ASSERT(frame.m_ptr);

	Reservoir & br = *(Reservoir*)(m_br + 0);

	FrameData & data = *(FrameData*)(m_framedata + 0);

	m_stream_ptr = frame.m_ptr;

	if (!Private::ReadSideInfo(*this, frame)) return 0;

	if (Private::ReadMain(*this, frame))	//returns false if not enough data in the bit reservoir
	{
		Private::DecodeFrame(*this, frame, data.is, output);

		return frame.m_length;		//constant 1152, if not first Info frame
	}
	else
	{
		return 0;
	}
}




//
//private impl

bool OpenMP3::Decoder::Private::ReadSideInfo(Decoder & self, const Frame & frame)
{
	FrameData & sideinfo = *(FrameData*)(self.m_framedata + 0);

	struct Local
	{
		FORCEINLINE static UInt32 GetSideBits(UInt n, const UInt8 * & ptr, UInt & idx)	//TODO unify with resorvoir ReadBits
		{
			UInt32 a;	// = *reinterpret_cast<const UInt32*>(ptr);

			UInt8 * p = reinterpret_cast<UInt8*>(&a);

			p[0] = ptr[3];
			p[1] = ptr[2];
			p[2] = ptr[1];
			p[3] = ptr[0];

			a = a << idx;		//Remove bits already used

			a = a >> (32 - n);	//Remove bits after the desired bits
									
			ptr += (idx + n) >> 3;

			idx = (idx + n) & 0x07;

			return a;
		}
	};

	bool mono = frame.m_mode == OpenMP3::kModeMono;

	UInt nch = (mono ? 1 : 2);

	UInt sideinfo_size = (mono ? 17 : 32);

	
	const UInt8 * ptr = ReadBytes(self, sideinfo_size);

	UInt idx = 0;


	sideinfo.main_data_begin = Local::GetSideBits(9, ptr, idx);

	Local::GetSideBits(mono ? 5 : 3, ptr, idx);	//skip private bits

												//scale factor selection information
	for (UInt ch = 0; ch < nch; ch++) for (UInt scfsi_band = 0; scfsi_band < 4; scfsi_band++) sideinfo.scfsi[ch][scfsi_band] = Local::GetSideBits(1, ptr, idx);

	for (UInt gr = 0; gr < 2; gr++)
	{
		for (UInt ch = 0; ch < nch; ch++)
		{
			sideinfo.part2_3_length[gr][ch] = Local::GetSideBits(12, ptr, idx);

			sideinfo.big_values[gr][ch] = Local::GetSideBits(9, ptr, idx);

			if (sideinfo.big_values[gr][ch] > 288) return false;

			sideinfo.global_gain[gr][ch] = Float32(Local::GetSideBits(8, ptr, idx));

			sideinfo.scalefac_compress[gr][ch] = Local::GetSideBits(4, ptr, idx);

			sideinfo.window_switching[gr][ch] = Local::GetSideBits(1, ptr, idx) == 1;

			if (sideinfo.window_switching[gr][ch])
			{
				sideinfo.block_type[gr][ch] = Local::GetSideBits(2, ptr, idx);

				sideinfo.mixed_block[gr][ch] = Local::GetSideBits(1, ptr, idx) == 1;

				for (UInt region = 0; region < 2; region++) sideinfo.table_select[gr][ch][region] = Local::GetSideBits(5, ptr, idx);

				for (UInt window = 0; window < 3; window++) sideinfo.subblock_gain[gr][ch][window] = Float32(Local::GetSideBits(3, ptr, idx));

				if ((sideinfo.block_type[gr][ch] == 2) && (!sideinfo.mixed_block[gr][ch]))
				{
					sideinfo.region0_count[gr][ch] = 8; /* Implicit */
				}
				else
				{
					sideinfo.region0_count[gr][ch] = 7; /* Implicit */
				}

				/* The standard is wrong on this!!! */   /* Implicit */
				sideinfo.region1_count[gr][ch] = 20 - sideinfo.region0_count[gr][ch];
			}
			else
			{
				for (UInt region = 0; region < 3; region++) sideinfo.table_select[gr][ch][region] = Local::GetSideBits(5, ptr, idx);

				sideinfo.region0_count[gr][ch] = Local::GetSideBits(4, ptr, idx);

				sideinfo.region1_count[gr][ch] = Local::GetSideBits(3, ptr, idx);

				sideinfo.block_type[gr][ch] = 0;  /* Implicit */
			}

			sideinfo.preflag[gr][ch] = Local::GetSideBits(1, ptr, idx);

			sideinfo.scalefac_scale[gr][ch] = Local::GetSideBits(1, ptr, idx);

			sideinfo.count1table_select[gr][ch] = Local::GetSideBits(1, ptr, idx);
		}
	}

	//UInt n = ptr - start;

	return true;
}

bool OpenMP3::Decoder::Private::ReadMain(Decoder & self, const Frame & frame)
{
	Reservoir & br = *(Reservoir*)(self.m_br + 0);

	FrameData & data = *(FrameData*)(self.m_framedata + 0);

	struct Local
	{
		static bool ReadMainData(Decoder & self, Reservoir & br, UInt main_data_begin, UInt main_data_size)
		{
			if (main_data_begin > br.main_data_top)
			{
				/* No,there is not,so we skip decoding this frame,but we have to
				* read the main_data bits from the bitstream in case they are needed
				* for decoding the next frame. */
				ReadBytes(self, main_data_size, &(br.main_data_vec[br.main_data_top]));

				/* Set up pointers */
				br.main_data_ptr = &(br.main_data_vec[0]);
				br.main_data_idx = 0;
				br.main_data_top += main_data_size;

				return false;    /* This frame cannot be decoded! */
			}

			/* Copy data from previous frames */
			for (UInt i = 0; i < main_data_begin; i++) br.main_data_vec[i] = br.main_data_vec[br.main_data_top - main_data_begin + i];

			/* Read the main_data from file */
			ReadBytes(self, main_data_size, br.main_data_vec + main_data_begin);

			/* Set up pointers */
			br.main_data_ptr = &(br.main_data_vec[0]);
			br.main_data_idx = 0;
			br.main_data_top = main_data_begin + main_data_size;

			return true;
		}
	};


	UInt nch = (frame.m_mode == OpenMP3::kModeMono ? 1 : 2);

	UInt sideinfo_size = (nch == 1 ? 17 : 32);

	UInt main_data_size = frame.m_datasize - sideinfo_size;


	//Assemble the main data buffer with data from this frame and the previous two frames into a local buffer used by the Get_Main_Bits function
	//main_data_begin indicates how many bytes from previous frames that should be used
	if (!Local::ReadMainData(self, br, data.main_data_begin, main_data_size)) return false; //This could be due to not enough data in reservoir

	br.hack_bits_read = 0;


	UInt sfb;

	for (UInt gr = 0; gr < 2; gr++)
	{
		auto scalefac_compress_gr = data.scalefac_compress[gr];

		auto scalefac_l_gr = data.scalefac_l[gr];

		auto scalefac_s_gr = data.scalefac_s[gr];

		for (UInt ch = 0; ch < nch; ch++)
		{
			UInt part_2_start = br.GetPosition();

			/* Number of bits in the bitstream for the bands */

			//UInt slen1 = kScaleFactorSizes[data.scalefac_compress[gr][ch]][0];

			//UInt slen2 = kScaleFactorSizes[data.scalefac_compress[gr][ch]][1];

			auto scalefactors = kScaleFactorSizes[scalefac_compress_gr[ch]];

			UInt slen1 = scalefactors[0];

			UInt slen2 = scalefactors[1];

			if (data.window_switching[gr][ch] && (data.block_type[gr][ch] == 2))
			{
				if (data.mixed_block[gr][ch]) for (sfb = 0; sfb < 8; sfb++) scalefac_l_gr[ch][sfb] = br.ReadBits(slen1);

				//for (sfb = 0; sfb < 12; sfb++)
				//{
				//	UInt nbits = (sfb < 6) ? slen1 : slen2; //TODO optimise, slen1 for band 3-5, slen2 for 6-11

				//	for (UInt win = 0; win < 3; win++) scalefac_s_gr[ch][sfb][win] = br.ReadBits(nbits);
				//}

				for (sfb = 0; sfb < 6; sfb++) for (UInt win = 0; win < 3; win++) scalefac_s_gr[ch][sfb][win] = br.ReadBits(slen1);

				for (sfb = 6; sfb < 12; sfb++) for (UInt win = 0; win < 3; win++) scalefac_s_gr[ch][sfb][win] = br.ReadBits(slen2);
			}
			else
			{ /* block_type == 0 if winswitch == 0 */
			  /* Scale factor bands 0-5 */
				if ((data.scfsi[ch][0] == 0) || (gr == 0))
				{
					for (sfb = 0; sfb < 6; sfb++) scalefac_l_gr[ch][sfb] = br.ReadBits(slen1);
				}
				else if ((data.scfsi[ch][0] == 1) && (gr == 1))
				{
					/* Copy scalefactors from granule 0 to granule 1 */
					for (sfb = 0; sfb < 6; sfb++) data.scalefac_l[1][ch][sfb] = data.scalefac_l[0][ch][sfb];
				}

				/* Scale factor bands 6-10 */
				if ((data.scfsi[ch][1] == 0) || (gr == 0))
				{
					for (sfb = 6; sfb < 11; sfb++) scalefac_l_gr[ch][sfb] = br.ReadBits(slen1);
				}
				else if ((data.scfsi[ch][1] == 1) && (gr == 1))
				{
					/* Copy scalefactors from granule 0 to granule 1 */
					for (sfb = 6; sfb < 11; sfb++) data.scalefac_l[1][ch][sfb] = data.scalefac_l[0][ch][sfb];
				}

				/* Scale factor bands 11-15 */
				if ((data.scfsi[ch][2] == 0) || (gr == 0))
				{
					for (sfb = 11; sfb < 16; sfb++) scalefac_l_gr[ch][sfb] = br.ReadBits(slen2);
				}
				else if ((data.scfsi[ch][2] == 1) && (gr == 1))
				{
					/* Copy scalefactors from granule 0 to granule 1 */
					for (sfb = 11; sfb < 16; sfb++) data.scalefac_l[1][ch][sfb] = data.scalefac_l[0][ch][sfb];
				}

				/* Scale factor bands 16-20 */
				if ((data.scfsi[ch][3] == 0) || (gr == 0))
				{
					for (sfb = 16; sfb < 21; sfb++) scalefac_l_gr[ch][sfb] = br.ReadBits(slen2);
				}
				else if ((data.scfsi[ch][3] == 1) && (gr == 1))
				{
					/* Copy scalefactors from granule 0 to granule 1 */
					for (sfb = 16; sfb < 21; sfb++) data.scalefac_l[1][ch][sfb] = data.scalefac_l[0][ch][sfb];
				}
			}

			data.count1[gr][ch] = ReadHuffman(br, frame.m_sr_index, data, part_2_start, gr, ch, data.is[gr][ch]);
		}
	}

	//TODO read ancillary data here

	//UInt bytes = (br.hack_bits_read / 8) + (br.hack_bits_read & 7 ? 1 : 0);

	//const UInt8 * ptr = frame.m_ptr + sideinfo_size + data.main_data_begin + bytes;

	//size_t size = (frame.m_ptr + frame.m_datasize) - ptr;

	return true;
}

void OpenMP3::Decoder::Private::DecodeFrame(Decoder & self, const Frame & header, Float32 datais[2][2][576], float out[2][1152])
{
	const FrameData & data = *(FrameData*)(self.m_framedata + 0);

	UInt sfreq = header.m_sr_index;

	auto & store = self.m_hs_store;

	auto & v_vec = self.m_sbs_v_vec;

	if (header.m_mode == OpenMP3::kModeMono)
	{
		for (UInt gr = 0; gr < 2; gr++)
		{
			auto & is = datais[gr][0];

			Requantize(sfreq, data, gr, 0, is);

			Reorder(sfreq, data, gr, 0, is);

			Antialias(data, gr, 0, is);

			HybridSynthesis(data, gr, 0, store[0], is);

			FrequencyInversion(is);

			SubbandSynthesis(data, is, v_vec[0], out[0] + (576 * gr));
		}

		memcpy(out[1], out[0], 1152 * sizeof(Float32));
	}
	else
	{
		UInt8 joint_stereo_mode = header.m_mode_extension;

		bool stereo_decoding = (header.m_mode == 1) && (joint_stereo_mode != 0);	//joint_stereo & (Intensity stereo | MS stereo)

		for (UInt gr = 0; gr < 2; gr++)
		{
			for (UInt ch = 0; ch < 2; ch++)
			{
				auto & is = datais[gr][ch];

				Requantize(sfreq, data, gr, ch, is);

				Reorder(sfreq, data, gr, ch, is);
			}

			if (stereo_decoding) Stereo(sfreq, joint_stereo_mode, data, gr, datais[gr]);

			for (UInt ch = 0; ch < 2; ch++)
			{
				auto & is = datais[gr][ch];

				Antialias(data, gr, ch, is);

				HybridSynthesis(data, gr, ch, store[ch], is);

				FrequencyInversion(is);

				SubbandSynthesis(data, is, v_vec[ch], out[ch] + (576 * gr));
			}
		}
	}
}

void OpenMP3::Decoder::Private::ReadBytes(Decoder & self, UInt no_of_bytes, UInt data_vec[])
{
	//TODO this should return pointer to bytes, not upscale to UInt32

	for (UInt i = 0; i < no_of_bytes; i++)
	{
		data_vec[i] = *self.m_stream_ptr++;
	}
}

FORCEINLINE const OpenMP3::UInt8 * OpenMP3::Decoder::Private::ReadBytes(Decoder & self, UInt nbytes)
{
	const UInt8 * ptr = self.m_stream_ptr;

	self.m_stream_ptr += nbytes;

	return ptr;
}