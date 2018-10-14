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

Program to classify input voices into one of the 5 vowels
It also generates cropped clips of the utterances based on a configurable threshold value(before compile)
and saves them after applying DC shift and Normalizing.

Input : input.txt containing the utterances
Output: A/E/I/O/U for each of utterances

*/

void writeLines(ofstream &OUTPUT_FS, vector<string> &lines) {
	for (int i = 0; i < lines.size(); ++i) {
		OUTPUT_FS << lines[i];
	}
}
void getLines(vector<string>& linesVec_, ifstream& inFile_)
{
	string line;
	while (getline(inFile_, line))
	{
		linesVec_.push_back(line);
	}
}

vector<double> hamming(vector<double> &samples,int s){
	vector<double> hamming_window(WINDOW_SIZE);
	for (int i = 0; i < WINDOW_SIZE; ++i)
		hamming_window[i]=samples[s+i]*(0.54 - 0.46 * cos( 2 * M_PI *  i / (WINDOW_SIZE - 1)));
	return hamming_window;
}

vector<double> getCis(vector<double> &Ais, double &E) {
	vector<double> Cis(Ais);
	Cis[0] = log2(E); 
	for (int m = 1; m <= P_ORDER; ++m) {
		for (double j = 1; j <= m - 1; ++j)
			Cis[m] += (j * Cis[j] * Ais[m - j])/(double)m;
	}
	return Cis;
}

vector<double> getAis(vector<double> &Ris) {
	vector<double> Ais(P_ORDER+1, 0);
	Ais[0] = 1.0;
	//ratio
	double K_iter, E_iter; 
	E_iter = Ris[0];
	for (int i = 1; i <= P_ORDER; ++i) {
		vector<double> prevAis(Ais);
		K_iter = Ris[i];
		for (int j = 1; j <= i - 1; ++j)K_iter -= (prevAis[j] * Ris[i - j]);
		// E_iter here is functionally prevE_iter
			K_iter /= E_iter;
		Ais[i] = K_iter;
		for (int j = 1; j <= i - 1; ++j) {
			Ais[j] = prevAis[j] - K_iter*prevAis[i - j];
		}
		E_iter *= (1 - K_iter*K_iter);
	}
	return Ais;
}
vector<double> getRis(vector<double> &window) {
	vector<double> Ris(P_ORDER + 1);
	for (int gap = 0; gap <= P_ORDER; ++gap) {
		for (int j = 0; j + gap < WINDOW_SIZE; ++j)
			Ris[gap] += window[j] * window[j + gap];
	}
	return Ris;
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

void saveCis(string filename, vector<vector<double> > &allCis, bool verbose=false) {
	int n=allCis.size(),p=allCis[0].size();
	ofstream OUTPUT_FS(filename, std::ios_base::app);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < p; ++j) 
			OUTPUT_FS << allCis[i][j] << "\t";
		OUTPUT_FS << endl;
	}
	OUTPUT_FS << endl;
	OUTPUT_FS.close();
	if (verbose)cout << "Appened to file: " << filename << endl;
}


#define FIXED_FLOAT(x) std::fixed <<std::setprecision(6)<<(x) 
int main(int argc, char const *argv[]) {
	char VOWELS[5] = { 'a','e','i','o','u' };
	for (int vowelIdx = 0; vowelIdx < 5; ++vowelIdx)
	{
		char TRAIN_VOWEL = VOWELS[vowelIdx];
		std::ostringstream filename;
		filename.str(std::string());
		filename << CURR_DIR<< INPUT_FOLDER <<TRAIN_VOWEL<<".txt";
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
			// record instead 
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
		if(INPUT_TRIMMED_ALREADY){
			clips.push_back(make_pair(0,N_SAMPLES));
			// clips.push_back(make_pair(N_SAMPLES/4,N_SAMPLES));
		}
		else{

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

		// For each clip perform the calculations
		filename.str(std::string());
		filename << "coeffs/cis_"<<TRAIN_VOWEL<<".txt";			
		//clear the file
		ofstream OUTPUT_FS(filename.str()); OUTPUT_FS.close();

		for (int i = 0; i < NUM_CLIPS; ++i) {
			clip_size = (clips[i].second - clips[i].first);
			cout << "Clip " << i << ": "
			<< "[" << FIXED_FLOAT(clips[i].first / SAMPLERATE) << "s," 
			<< FIXED_FLOAT(clips[i].second / SAMPLERATE) << "s] "
			<< "\tDuration: " << FIXED_FLOAT(clip_size / SAMPLERATE) << "s (" << clip_size << " samples)\n";
			
			vector<vector<double> > allCis;
			// better region check run windows on from offset of clip_size/4
			for (int s = clips[i].first+clip_size/4, w = 1; s + WINDOW_SIZE < clips[i].second; s += WINDOW_STRIDE, w++) {
				vector<double> windowCopy;
				if(USE_HAMMING)
				//apply hamming window at s			
					windowCopy = hamming(samples,s);
				else
				// without hamming window
					windowCopy.insert(windowCopy.end(),samples.begin()+s,samples.begin()+s+WINDOW_SIZE);

				cout << "Window " << w << " [" << s << "," <<(s+WINDOW_SIZE)<<"] :\t";
				vector<double> Ris = getRis(windowCopy);
				if(SHOW_RIs){
					cout << "Ris: ";
					for (int i = 0; i <= P_ORDER; ++i)
						cout << Ris[i] << " ";
					cout << endl;
				}
				vector<double> Ais = getAis(Ris);
				if(SHOW_AIs){
					cout << "Ais: ";
					for (int i = 0; i <= P_ORDER; ++i)
						cout << Ais[i] << " ";
					cout << endl;
				}
				vector<double> Cis = getCis(Ais, Ris[0]);
				if(SHOW_CIs){
					cout << "Cis: ";
					for (int i = 0; i <= P_ORDER; ++i)
						cout << Cis[i] << " ";
					cout << endl;
				}
				// cout << endl;
				allCis.push_back(Cis);
				if(w==FIRST_N_FRAMES) break;
			}

			saveCis(filename.str(), allCis,true);
			cout << endl;
		}
		if(SAVE_TRIMMED_CLIPS){
			cout<<"Generating filenames(for saving clips) for vowel: "<<TRAIN_VOWEL<<endl;
			for (int i = 0; i < NUM_CLIPS; ++i) {
				// clear the stream;
				filename.str(std::string());
				filename << "output_clips/150101021_"<<TRAIN_VOWEL<<"_"<<i<<".txt";			
				//write cropped sample into new file		
				writeSample(filename.str(),samples,clips[i].first,clips[i].second);
			}
			cout<<"All clips saved in 'output_clips/' directory\n";
		}
	}

	return 0;
}