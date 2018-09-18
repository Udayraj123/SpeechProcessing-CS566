## Description
This is a program to classify input voices as One or Six
It also generates cropped clips of the utterances, based on a threshold value(configurable via a variable).

## How to run the code
1. Compile main.cpp (Important: use the '-std=c++11' option else it wont compile)
	g++ -std=c++11 -o main main.cpp
2. Run main.exe
	main.exe  <duration_in_seconds>
3. When it shows 'Start Recording...' you can start speaking. You can speak many utterances until it shows 'Stop Recording.'. Then press any key to continue. The program will separate them and tell for each of them if it was a ONE or a SIX.
All the utterances found will be saved in the 'output_clips/' folder with their labels.

## Providing your input
1. By default the program will take your Voice as the input for given duration
	main.exe 5
The recording is also saved in the file input.txt for later use.

2. If no argument is provided, the program reads from 'input.txt' instead of recording. 
	main.exe
where input.txt is present in current folder
