*** OpenMP3 v0.95 ***

- A licence-free MP3 decoding library.




* CREDITS & REFERENCES * 

- PDMP3 is a licence-free MP3 decoding library written by Krister Lagerström (krister@kmlager.com)

- OpenMP3, by AudioBoy, is essentially a refactor / cleanup of this library

- See the original PDMP3 library here: https://github.com/technosaurus/PDMP3

- "PDMP3-5d42a7951d977ab8a6f54657ed96e49cf532a946/" is a snapshot of PDMP3
	- this snapshot was used for the basis of OpenMP3 
	- it is intentionally not the latest commit of PDMP3, made by adalinbv, which adding a streaming API as that latest commit had buffer overrun issues which were unfixed at the time of writing OpenMP3




* LIBRARY DESIGN CONSIDERATIONS * 

- Clean, modern C++ interface

- No external dependencies (not even to C or C++ standard includes)

- No memory allocations inside the library

- Multi-instance / thread-safe (main technical limitation of PDMP3)

- Proper seperation of concerns (decoupled frame streaming & parsing from the decoding algorithm)

- Single .h and .cpp include




* TODO * 
	
- Support for iterating / extracting ID3 tags
	- Currently they should just be skipped

- Optimisations
	- Some basic optimisations were done as part of the refactor
	- However lots of optimisation potential still exists.  Roughly three blocks:
		- (1) re-factor data layout to reduce deferencing
		- (2) some more calculations to be pre-computed (maybe)
		- (3) use simd/vector operations




* CHANGE LOG * 

- Version 0.95
	- Fixed an error in kHuffmanTables, seems to have fixed known decoding errors
	
- Version 0.91
	- Approx 3x faster decoding due to use of pre-calc tables in IMDCT_Win
	- Skip initial silent 'Info' frame created by LAME encoder
	
- Version 0.9: Initial release

