
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
double tokhuraDist(vector <double> cis,vector <double> cis2){
	static double weights[]={1,3,7,13,19,22,25,33,42,50,56,61};
	double d = 0 ;
	for (int i = 1; i < cis.size(); ++i){		
		d += (cis[i]-cis2[i])*(cis[i]-cis2[i])*weights[i-1];
	}
	return d;
}
double vectorDist(vector <double> cis,vector <double> cis2){
	double d = 0 ;
	for (int i = 1; i < cis.size(); ++i){
		d += (cis[i]-cis2[i])*(cis[i]-cis2[i]);
	}
	return d;
}
double calcAvgDistance(vector <double> cis,vector < vector <double> >frameCis ){
	double avg = 0 , d;
	for (int f = 0; f < frameCis.size(); ++f) {
		if(USE_TOKHURA)
			d = tokhuraDist(cis,frameCis[f]);
		else
			d = vectorDist(cis,frameCis[f]);
		avg += d;
	}
	avg /= frameCis.size();
	return avg;
}
vector <vector< vector<double> > >  loadFrameWiseCis(string filename){
	ifstream INPUT_FS;
	INPUT_FS.open(filename);
	vector <vector< vector<double> > > frameWiseCis(FIRST_N_FRAMES);
	int frame=0;
	// cout<<filename<<endl;
	while(true){
		vector<double>temp(P_ORDER+1);
		for (int p = 0; p <= P_ORDER; ++p){
			INPUT_FS>>temp[p];			
		}
		if(INPUT_FS.eof())break;
		frameWiseCis[frame].push_back(temp);
		frame = (frame+1)%FIRST_N_FRAMES;		
	}
	// cout<<filename<<endl;
	INPUT_FS.close();
	return frameWiseCis;
}

#define FIXED_FLOAT(x) std::fixed <<std::setprecision(6)<<(x) 
int main(int argc, char const *argv[]) {
	std::ostringstream filename;
	filename.str(std::string());
	if(argc>1){
		filename << argv[1];
	}
	else{
		filename << CURR_DIR <<"_";
	}
	ifstream INPUT_FS;
	INPUT_FS.open(filename.str());

	if(!INPUT_FS.good()) { 
		// .good() returns if file not exists
		string duration = DEFAULT_DURN;
		filename.str(std::string());
		filename << CURR_DIR<< DEFAULT_RECORD_FILE;
		cout << "No/Invalid input file provided, Calling recording module\n Duration:  " << duration << "s" << endl;
		cout <<" Note: You can speak more than one utterances\n";
		system(("Recording_Module.exe " + duration + " temp.wav " + filename.str()).c_str());
		// calling open on an already open stream fails//clear flags too
		INPUT_FS.close(); INPUT_FS.clear();
		INPUT_FS.open(filename.str());
		remove("temp.wav");
		// cout<<"ERROR: Invalid Input file: "<<filename.str()<<endl;
		// return 1;
	}
	char VOWELS[5] = { 'a','e','i','o','u' };
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
	if(TEST_TRIMMED_ALREADY){
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
			}
				// else default is false
		}
			// Detect clips using applied threshold
		bool prev = blurred_samples[0];
		int start = 0;
		for (int i = 1; i < N_SAMPLES; ++i) {
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
	
	// read stored Ci values for each vowel
	vector< vector< vector<double> > >	vowelWiseFrameCis [5];
	for (int vowelIdx = 0; vowelIdx < 5; ++vowelIdx)
	{
		filename.str(std::string());
		filename << CURR_DIR << "coeffs/cis_"<<VOWELS[vowelIdx]<<".txt";			
		vowelWiseFrameCis[vowelIdx] = loadFrameWiseCis(filename.str());
	}
	
	// For each clip perform the calculations
	for (int i = 0; i < NUM_CLIPS; ++i) {
		clip_size = (clips[i].second - clips[i].first);
		cout << "Clip " << i << ": "
		<< "[" << FIXED_FLOAT(clips[i].first / SAMPLERATE) << "s," 
		<< FIXED_FLOAT(clips[i].second / SAMPLERATE) << "s] "
		<< "\tDuration: " << FIXED_FLOAT(clip_size / SAMPLERATE) << "s (" << clip_size << " samples)\n";

		vector<vector<double> > frameWiseCis;
		// better region check run windows on from offset of clip_size/4
		for (int s = clips[i].first + clip_size/4, w = 1; s + WINDOW_SIZE < clips[i].second; s += WINDOW_STRIDE, w++) {
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
			frameWiseCis.push_back(Cis);
			if(w==FIRST_N_FRAMES) break;
		}

		/// Print distances
		double minAvg=100000;
		char foundVowel='X';
		cout << "Vowel \t Distance\t\n";
		for (int vowelIdx = 0; vowelIdx < 5; ++vowelIdx)
		{
			vector<double> vowelDistances(FIRST_N_FRAMES,0);
			double avgDistance=0;
			for (int w = 0; w < 5; ++w){
				vowelDistances[w] = calcAvgDistance(frameWiseCis[w],vowelWiseFrameCis[vowelIdx][w]);
				avgDistance += vowelDistances[w];
			}
			// average it
			avgDistance /= FIRST_N_FRAMES;
			if(minAvg > avgDistance){
				minAvg = avgDistance;
				foundVowel=VOWELS[vowelIdx];
			}
			// for (int w = 0; w < FIRST_N_FRAMES; ++w){
			// 	cout<<VOWELS[vowelIdx]<<"\t F"<<w<<"\t "<<vowelDistances[w]<<"\n";
			// }
			cout<<VOWELS[vowelIdx]<<"\t "<<avgDistance<<endl;
		}
		cout << "\n\tDetected Vowel : \'"<<foundVowel<<"\' at minimum avg distance of "<<minAvg<<endl<<endl;
	}

	return 0;
}