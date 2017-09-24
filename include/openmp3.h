#pragma once




//
//declarations

namespace OpenMP3
{

	//enums
	
	enum Mode
	{
		kModeStereo,
		kModeJointStereo,
		kModeDualMono,
		kModeMono
	};



	//classes

	class Library;

	class Iterator;					//iterate an in memory mp3 file

	//TODO class StreamIterator;	//has internal buffer, used for real-time streaming

	class Decoder;

	class Frame;



	//integral types

	typedef unsigned char UInt8;

	typedef unsigned short UInt16;

	typedef unsigned int UInt32;

	typedef UInt32 UInt;


	typedef float Float32;

};





//
//library

class OpenMP3::Library
{
public:

	Library();



private:

	friend Decoder;

};





//
//iterator

class OpenMP3::Iterator
{
public:

	//lifetime

	Iterator(const Library & library, const UInt8 * data, UInt size);



	//access

	bool GetNext(Frame & frame);



private:

	struct Private;


	const UInt8 * m_start, * m_ptr, * m_end;

	bool m_hack_first;		//for Info frame skipping hack

};




//
//decoder

class OpenMP3::Decoder
{
public:

	//lifetime
	
	Decoder(const Library & library);

	
	
	//access
	
	void Reset();															//Call before processing new song, resets internal buffers

	UInt ProcessFrame(const Frame & frame, Float32 output[2][1152]);		//return: number samples extracted from frame (0 or 1152)



private:

	struct Private;


	const Library & library;


	UInt8 m_br[8320];

	UInt8 m_framedata[10496];

	Float32 m_hs_store[2][32][18];

	Float32 m_sbs_v_vec[2][1024];


	const UInt8 * m_stream_ptr;

	UInt m_stream_remainder;

};




//
//frame

class OpenMP3::Frame
{
public:

	//lifetime

	Frame();



	//access

	//TODO UInt GetPosition() const;

	UInt GetBitRate() const;

	UInt GetSampleRate() const;

	Mode GetMode() const;



private:

	friend Iterator;

	friend Decoder;


	UInt8 m_bitrate_index;

	UInt8 m_sr_index;

	Mode m_mode;

	UInt8 m_mode_extension;


	const UInt8 * m_ptr;		//pointer to data area

	UInt m_datasize;			//size of whole frame, minus headerword + check

	UInt m_length;				//for Info frame skipping

};
