#include "config.h"
#include <iostream>
#include <string>
#include <vector>
#include <set> // debug fn names list
#include <algorithm> //reverse and sort
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

#define FIXED_FLOAT(x) std::fixed <<std::setprecision(6)<<(x) 
set <string> DEBUG_FNs;
#define DEBUG(x) if(DEBUG_ON && DEBUG_FNs.count(__FUNCTION__)){cout << ">> " << #x << " : \t";cout<<FIXED_FLOAT(x)<<endl;}
#define DEBUGV(x) if(DEBUG_ON && DEBUG_FNs.count(__FUNCTION__)){cout << ">> " << #x << " : \t";for(auto _i = (x).begin();_i!=(x).end();_i++)cout<<FIXED_FLOAT(*_i)<<"\t";cout<<"\n";}
#define PRINT(x) {cout << #x << " : \t";cout<<FIXED_FLOAT(x)<<endl;}
#define PRINTV(x){cout << #x << " : \t";for(auto _i = (x).begin();_i!=(x).end();_i++)cout<<FIXED_FLOAT(*_i)<<" ";cout<<"\n";}
#define FOR(i,n)for (ll i = 0; i < (n); ++i)
#define FORG(i,a,b)for (ll i = (a); i <= (b); ++i)
#define FORD(i,a,b)for (ll i = (a); i >= (b); --i)

#define pb push_back
#define vi vector <int> 
#define vvi vector <vector<int> >
#define vd vector <double> 
#define vvd vector <vector<double> >
#define vvvd vector <vector<vector<double> > >

/*

Program to generate universe file from the trimmed recordings of 10 digits

Input : digitdata folder containing the digitwise trimmed utterances 
Output: The universe file as well as individual files having Ci values of each digit

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

void saveCis(string filename, vector<vector<double> > &allCis, bool verbose=false) {
	int n=allCis.size(),p=allCis[0].size();
	ofstream OUTPUT_FS(filename, std::ios_base::app);
	for (int i = 0; i < n; ++i) {
		// skip first
		// for (int j = 1; j < p; ++j) 
		for (int j = 0; j < p; ++j) 
			OUTPUT_FS << allCis[i][j] << "\t";
		OUTPUT_FS << endl;
	}
	OUTPUT_FS << endl;
	OUTPUT_FS.close();
	if (verbose)cout << "Appended to file: " << filename << endl;
}

int main(int argc, char const *argv[]) {
	DEBUG_FNs.insert("main");
	char DIGITS[] = { '0', '1','2','3','4','5','6','7','8','9' };
	std::ostringstream filename;
	//clear the stream
	filename.str(std::string());
	filename << CURR_DIR << "coeffs/universe.txt";
	string UNIVERSE_FILE = filename.str();
	//clear the file
	ofstream OUTPUT_FS(filename.str()); OUTPUT_FS.close();
	FOR(UTTERNO,10){
		for (int digitIdx = 0; digitIdx < sizeof(DIGITS)/sizeof(DIGITS[0]); ++digitIdx)
		{
			char TRAIN_DIGIT = DIGITS[digitIdx];
			filename.str(std::string());
			filename << CURR_DIR<< "digit_data1/150101021_"<<TRAIN_DIGIT<<"_"<<UTTERNO<<".txt";
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
			clips.push_back(make_pair(0,N_SAMPLES));
		/*
			Now at this point input is dc shifted(if any), normalized, and trimmed into different clips
		*/
			int clip_size, NUM_CLIPS = clips.size();
			cout << NUM_CLIPS << " clips found\n";
			if (NUM_CLIPS > N_CLIPS_LIMIT) {
				cout << "Error: too many clips found: " << NUM_CLIPS << endl;			
				return 0;
			}
			else if (NUM_CLIPS == 0) {
				cout << "Error: No clips found: either increase recording volume or check if mic is working.\n";
				return 0;
			}
			cout << "Total duration: " << FIXED_FLOAT(N_SAMPLES / SAMPLERATE) << "s (" << N_SAMPLES << " samples)" << endl;
			
			//clear the stream
			filename.str(std::string());
			filename << CURR_DIR << "coeffs/"<<TRAIN_DIGIT<<"_"<<UTTERNO<<"_cis.txt";
			// clear the individual cis file
			ofstream OUTPUT_FS(filename.str()); OUTPUT_FS.close();

			for (int i = 0; i < NUM_CLIPS; ++i) {
				clip_size = (clips[i].second - clips[i].first);
				cout << "Clip " << i << ": "
				<< "[" << FIXED_FLOAT(clips[i].first / SAMPLERATE) << "s," 
				<< FIXED_FLOAT(clips[i].second / SAMPLERATE) << "s] "
				<< "\tDuration: " << FIXED_FLOAT(clip_size / SAMPLERATE) << "s (" << clip_size << " samples)\n";

				vector<vector<double> > allCis;
			// better region check run windows on from offset of clip_size/4
				// TODO: tune this factor
				int CLIP_MARGIN = clip_size/4;
				for (int s = clips[i].first + CLIP_MARGIN, w = 1; s + WINDOW_SIZE < clips[i].second - CLIP_MARGIN; s += WINDOW_STRIDE, w++) {
					
					vector<double> windowCopy;
					// whether to apply a hamming window at s			
					if(USE_HAMMING)
						windowCopy = hamming(samples,s);
					else
						windowCopy.insert(windowCopy.end(),samples.begin()+s,samples.begin()+s+WINDOW_SIZE);

					cout << "Window " << w << " [" << s << "," <<(s+WINDOW_SIZE)<<"] :\t";
					
					vector<double> Ris = getRis(windowCopy);
					if(SHOW_RIs)
						PRINTV(Ris)

					vector<double> Ais = getAis(Ris);
					if(SHOW_AIs)
						PRINTV(Ais)

					vector<double> Cis = getCis(Ais, Ris[0]);
					if(SHOW_CIs)
						PRINTV(Cis)

					allCis.push_back(Cis);

					if(w==OBSERVATIONS_LIM) 
						break;
				}

				saveCis(UNIVERSE_FILE, allCis,true);
				saveCis(filename.str(), allCis,true);
				cout << endl;
			}
		}
	}

	return 0;
}