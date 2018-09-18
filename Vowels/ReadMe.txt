## Description
This is a program to classify input voices into one of the 5 vowels
It also generates cropped clips of the utterances based on a configurable threshold value(before compile) and saves them after applying DC shift and Normalizing.

## How to run the code
0. Set config.h file: Please provide following inputs before proceeding to compile the program
	CURR_DIR - Full path to current directory of this program (Including the '/' at the end)
	
	INPUT_FOLDER - If you want to regenerate Ci files for your recordings. Put them in this directory. It should have five files present: a.txt, e.ext, i.txt.o.txt and u.txt
	INPUT_TRIMMED_ALREADY - whether the above five files are already trimmed
	TEST_TRIMMED_ALREADY - whether the test file is already trimmed
	
	SHOW_RIs - Whether to print Ri values on terminal
	SHOW_AIs - Whether to print Ai values on terminal
	SHOW_CIs - Whether to print Ci values on terminal
	
	Model parameters:
	USE_HAMMING 	 1 for using hamming, 0 for using step
	USE_TOKHURA	 	 1 for using tokhura, 0 for using euclidean
	Further options are provided in the config.h file

1. Compile recognize_vowels.cpp (Important: use the '-std=c++11' option else it wont compile)
	g++ -std=c++11 -o recognize_vowels recognize_vowels.cpp

2. Run recognize_vowels.exe
	recognize_vowels.exe  <filename>
If <filename> is the test file to recognize vowel utterances from.
(If the file is trimmed already, please change TEST_TRIMMED_ALREADY to 1 in config.h before compiling)
After the file is processed, the program will output the Tokhura's(or euclidean depending on config above) distances for each utterance in the input
And the one with minimum distance is shown as detected vowel.

If <filename> above is not present, the recording will start with duration DEFAULT_DURN (set from config.h file).
When it shows 'Start Recording...' you can start speaking. You can speak many utterances until it shows 'Stop Recording.'. Then press any key to continue.
