#pragma once

#include "../include/openmp3.h"




//
//external dependencies

#include <math.h>
#include <cstring>	//for memset, memcopy




//
//macros

#define CT_ASSERT(value) {int _ct_assert[value] = {0};}

#ifndef ASSERT
#define ASSERT(x)
#endif 

#ifndef FORCEINLINE
#ifdef __APPLE__
#define FORCEINLINE __attribute((always_inline)) inline
#elif _WIN32
#define FORCEINLINE __forceinline
#elif _WIN64
#define FORCEINLINE __forceinline
#else
#define FORCEINLINE inline
#endif
#endif




//
//declarations

namespace OpenMP3
{

	enum Layer
	{
		kLayerReserved,
		kLayer3,
		kLayer2,
		kLayer1
	};


	struct Reservoir;

	struct FrameData;


	typedef double Float64;
	


	//helper functions

	inline void MemClear(void * address, size_t size)
	{
		memset(address, 0, size);
	}

}




//
//reservoir

struct OpenMP3::Reservoir
{
	void SetPosition(UInt bit_pos);

	UInt GetPosition();

	UInt ReadBit();

	UInt ReadBits(UInt n);

	UInt main_data_vec[2 * 1024];	//Large data	TODO bytes
	UInt *main_data_ptr;			//Pointer into the reservoir
	UInt main_data_idx;				//Index into the current byte(0-7)
	UInt main_data_top;				//Number of bytes in reservoir(0-1024)

	UInt hack_bits_read;
};




//
//framedata

struct OpenMP3::FrameData
{
	unsigned main_data_begin;         /* 9 bits */


	//side

	unsigned scfsi[2][4];             /* 1 bit */
	unsigned part2_3_length[2][2];    /* 12 bits */
	unsigned big_values[2][2];        /* 9 bits */
	Float32 global_gain[2][2];       /* 8 bits */
	unsigned scalefac_compress[2][2]; /* 4 bits */

	bool window_switching[2][2];
	bool mixed_block[2][2];

	unsigned block_type[2][2];        /* 2 bits */
	unsigned table_select[2][2][3];   /* 5 bits */
	Float32 subblock_gain[2][2][3];  /* 3 bits */
									  /* table_select[][][] */
	unsigned region0_count[2][2];     /* 4 bits */
	unsigned region1_count[2][2];     /* 3 bits */
									  /* end */
	unsigned preflag[2][2];           /* 1 bit */
	unsigned scalefac_scale[2][2];    /* 1 bit */
	unsigned count1table_select[2][2];/* 1 bit */

	unsigned count1[2][2];            //calc: by huff.dec.!


	//main

	unsigned scalefac_l[2][2][21];    /* 0-4 bits */
	unsigned scalefac_s[2][2][12][3]; /* 0-4 bits */

	float is[2][2][576];               //calc: freq lines
};




//
//impl

inline void OpenMP3::Reservoir::SetPosition(UInt bit_pos)
{
	main_data_ptr = &(main_data_vec[bit_pos >> 3]);

	main_data_idx = bit_pos & 0x7;
}

FORCEINLINE OpenMP3::UInt OpenMP3::Reservoir::GetPosition()
{
	UInt pos = ((size_t)main_data_ptr) - ((size_t) &(main_data_vec[0]));

	return (pos * 2) + main_data_idx;
}

FORCEINLINE OpenMP3::UInt OpenMP3::Reservoir::ReadBit()
{
	hack_bits_read++;

	UInt tmp = main_data_ptr[0] >> (7 - main_data_idx);

	tmp &= 0x01;

	main_data_ptr += (main_data_idx + 1) >> 3;

	main_data_idx = (main_data_idx + 1) & 0x07;

	return tmp;
}

FORCEINLINE OpenMP3::UInt OpenMP3::Reservoir::ReadBits(UInt n)	//number_of_bits to read(max 24).  bits are returned in the LSB of the return value
{
	hack_bits_read += n;

	if (n == 0) return 0;

	/* Form a word of the next four bytes */
	UInt tmp = (main_data_ptr[0] << 24) | (main_data_ptr[1] << 16) | (main_data_ptr[2] << 8) | (main_data_ptr[3] << 0);

	/* Remove bits already used */
	tmp = tmp << main_data_idx;

	/* Remove bits after the desired bits */
	tmp = tmp >> (32 - n);

	/* Update pointers */
	main_data_ptr += (main_data_idx + n) >> 3;
	main_data_idx = (main_data_idx + n) & 0x07;

	return tmp;
}
