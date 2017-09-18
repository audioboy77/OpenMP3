*** OpenMP3 v0.9 ***

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

- Proper seperation of concerns (primarily, fully decoupling streaming & parsing from the actual decoding algorithm)




* KNOWN ISSUES * 

- Decoding errors exist!!! (as with original PDMP3)
	- See "fail/" folder for examples
	- I am *not* going to fix, I am not a DSP developer! Contributors are needed to fix this!


	

* TODO * 
	
- Optimisations
	- Some basic optimisations were done as part of the refactor
	- However lots of optimisation potential still exists.  Roughly two major blocks:
	- Phase (1) re-factor data layout to massively reduce deferencing, move some local static data to the Library class, plus other small optimisations 
	- Phase (2) some calculations to be pre-computed
	- Phase (3) use vector operations

