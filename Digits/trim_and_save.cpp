#include "config.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> // for i in filename
#include <iomanip> // for setprecision
# define M_PI           3.14159265358979323846
#include <math.h> //for log2,abs()
#include <stdlib.h> // for system()
using namespace std;

#define ll long long int
#ifdef DEBUG_MODE
#define DEBUG_ON true
#else
#define DEBUG_ON false
#endif
#define DEBUGV(x) if(DEBUG_ON){cout << ">> " << #x << " : \t";for(auto i = (x).begin();i!=(x).end();i++)cout<<(*i)<<"\t";cout<<"\n";}

/*

Trimming out the utterances

*/

void getLines(vector<string>& linesVec_, ifstream& inFile_)
{
	string line;
	while (getline(inFile_, line))
	{
		linesVec_.push_back(line);
	}
}
void writeSample(string filename, vector<double> &samples, int &s, int &e, bool verbose = false) {
	ofstream OUTPUT_FS(filename);
	OUTPUT_FS << "SAMPLES: " << (e - s + 1) << endl;
	OUTPUT_FS << "BITSPERSAMPLE:	16" << endl;
	OUTPUT_FS << "CHANNELS:	1" << endl;
	OUTPUT_FS << "SAMPLERATE:	" << SAMPLERATE << endl;
	OUTPUT_FS << "NORMALIZED:	FALSE" << endl;
	for (int i = s; i <= e; ++i) {
		OUTPUT_FS << samples[i] << endl;
	}
	OUTPUT_FS.close();
	if (verbose)cout << "Saved file: " << filename << endl;
}
#define FIXED_FLOAT(x) std::fixed <<std::setprecision(6)<<(x) 
int main(int argc, char const *argv[]) {
	// char DIGITS[5] = { 'a','e','i','o','u' };
	char DIGITS[] = { '0','1','2','3','4','5','6','7','8','9' };
	std::ostringstream filename;
	//clear the stream
	filename.str(std::string());
	for (int digitIdx = 0; digitIdx < sizeof(DIGITS)/sizeof(DIGITS[0]); ++digitIdx)
	{
		char TRAIN_DIGIT = DIGITS[digitIdx];
		filename.str(std::string());
		// filename << CURR_DIR<< INPUT_FOLDER <<TRAIN_DIGIT<<"_1.txt";
		// filename << CURR_DIR<< INPUT_FOLDER <<TRAIN_DIGIT<<"_2.txt";
		filename << CURR_DIR<< INPUT_FOLDER <<TRAIN_DIGIT<<"_3.txt";
		ifstream INPUT_FS;
		INPUT_FS.open(filename.str());

		if(!INPUT_FS.good()) { 
		// .good() returns if file not exists
			cout<<"ERROR: Invalid Input file: "<<filename.str()<<endl;
			return 1;
		}
	//take file input
		vector <string> lines;
		getLines(lines, INPUT_FS);
	//file loaded in memory, can close it now.
		INPUT_FS.close();

		if (!lines.size()) {
			cout << "ERROR: Empty data! Please check if correct filename is provided and try again.\n";
			return 1;
		}
		// skip first text lines 
		int OFFSET=0;
		while (!(lines[OFFSET][0]=='-' || lines[OFFSET][0]<='9' && lines[OFFSET][0]>='0')) {
			OFFSET++;
		}
		if(OFFSET>5){
			cout<<"Warning: File seems phishy!"<<endl;
		}

		// Normalization variables
		int N_SAMPLES = lines.size() - OFFSET;
		vector<double> samples;
		double x,maxX=0;

		//Read line by line, do DC shift as well as load the samples
		for (int i = 0; i < N_SAMPLES; ++i) {
			sscanf(&lines[i + OFFSET][0], "%lf", &x);
			x -= DC_SHIFT;
			samples.push_back(x);
			maxX = maxX<(double)abs(x) ? (double)abs(x) : maxX;
		}
		//Normalize
		double N_fac = N_AMP / maxX;

		if (DEBUG_ON)
			cout<<N_fac<<" "<<maxX<<endl;

		for (int i = 0; i < N_SAMPLES; ++i)
			samples[i] *= N_fac; 

		vector< pair<int, int> > clips;
		// do Trimming
		{
			vector<bool> blurred_samples(N_SAMPLES,false);
			double windowSum = 0;
			if (DEBUG_ON)
				cout<<windowSum<<" "<<M_SILENCE_ENERGY<<endl;
			// Apply M_SILENCE_ENERGY for windowed threshold on square of amplitude
			for (int i = 0; i < N_SAMPLES; ++i) {
				//x = abs(samples[i]);
				x = pow(samples[i], 2);
				windowSum += x;
				if (i>SILENCE_WINDOW) {
					//abs(samples[i-SILENCE_WINDOW]);
					windowSum -= pow(samples[i - SILENCE_WINDOW], 2);
					if ((windowSum / (double)SILENCE_WINDOW)>M_SILENCE_ENERGY)
						blurred_samples[i - SILENCE_WINDOW / 6] = true;
					// need another blurred_samples for correcting right end limit [think]
				}
				// else default is false
			}
			// Detect clips using applied threshold
			bool prev = blurred_samples[0];
			int start = 0; // offset from start
			for (int i = start+1; i < N_SAMPLES; ++i) {
				if (prev != blurred_samples[i]) {
					if (!start) {
						start = i;
					}
					else {
				//end
						clips.push_back(make_pair(start, i));
						start = 0;
					}
				}
				prev = blurred_samples[i];
			}
		}
		/*
			Now at this point input is dc shifted(if any), normalized, and trimmed into different clips
		*/
		int clip_size, NUM_CLIPS = clips.size();
		cout << NUM_CLIPS << " clips found\n";
		if (NUM_CLIPS > N_CLIPS_LIMIT) {
			// string filen = "blurred_samples.txt";
			// int s = 0;
			// writeSample(filen, blurred_samples, s, N_SAMPLES);
			cout << "Error: too many clips found: " << NUM_CLIPS << endl;
			cout << "Either increase limit or check blurred_samples.txt to resolve the issue\n";
			return 0;
		}
		else if (NUM_CLIPS == 0) {
			cout << "Error: No clips found: either increase recording volume or check if mic is working.\n";
			return 0;
		}
		cout << "Total duration: " << FIXED_FLOAT(N_SAMPLES / SAMPLERATE) << "s (" << N_SAMPLES << " samples)" << endl;
		
		cout<<"Generating filenames(for saving clips) for digit: "<<TRAIN_DIGIT<<endl;
		for (int i = 0; i < NUM_CLIPS; ++i) {
				// clear the stream;
			filename.str(std::string());
			filename << CURR_DIR << "output_clips/150101021_"<<TRAIN_DIGIT<<"_"<<i<<".txt";			
				//write cropped sample into new file		
			writeSample(filename.str(),samples,clips[i].first,clips[i].second);
		}
		cout<<"All clips saved in 'output_clips/' directory\n";
	}

	return 0;
}