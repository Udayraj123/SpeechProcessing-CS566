## How to run the code
1. Compile HMM.cpp 
	Important: use the '-std=c++11' option else it wont compile
	g++ -std=c++11 -o HMM.exe HMM.cpp

2. Run HMM.exe
	HMM.exe  

It will first output the P(O/lambda) according to forward procedure.

It will then output the PStar,QStar and Expected Observation Sequence according to Viterbi algorithm.