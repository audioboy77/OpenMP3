#include "../openmp3.h"
#include <vector>




//--- EXAMPLE #1 ---
//parsing through an mp3 file, to test if all frames are mono and get max bitrate

void ParseMP3(const void * mp3_file, unsigned file_size)
{
	OpenMP3::Library openmp3;
	
	OpenMP3::Iterator itr(openmp3, (OpenMP3::UInt8*)mp3_file, file_size);


	bool mono = true;

	OpenMP3::UInt br = 0;
	

	OpenMP3::Frame frame;

	while (itr.GetNext(frame))
	{
		if (frame.GetBitRate() > br) br = frame.GetBitRate();

		mono = mono && (frame.GetMode() == OpenMP3::kModeMono);
	}
}




//--- EXAMPLE #2 ---
//loading an mp3 file
//for illustration only! appending samples one-by-one to a std::vector is ineffcient...

void LoadMP3(const void * mp3_file, unsigned file_size)
{
	OpenMP3::Library openmp3;

	OpenMP3::Iterator itr(openmp3, (OpenMP3::UInt8*)mp3_file, file_size);


	float buffer[2][1152];

	std::vector <float> channels[2];

	bool mono = true;


	OpenMP3::Decoder decoder(openmp3);

	OpenMP3::Frame frame;

	while (itr.GetNext(frame))
	{
		OpenMP3::UInt nsamples = decoder.ProcessFrame(frame, buffer);

		for(int ch = 0; ch < 2; ++ch)
		{
			auto & channel = channels[ch];

			auto * data = buffer[ch];

			for (OpenMP3::UInt idx = 0; idx < nsamples; ++idx) channel.push_back(*data++);
		}

		mono = mono && (frame.GetMode() == OpenMP3::kModeMono);
	}

	OpenMP3::UInt length = channels[0].size();
}
