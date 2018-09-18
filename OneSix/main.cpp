#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> // for i in filename
#include <iomanip> // for setprecision
#include <cmath> //for abs()
#include <stdlib.h> // for system()
using namespace std;


/*

Program to classify input voices into one or six
It also generates cropped clips of the utterances based on a configurable threshold value(before compile) 
and saves them after applying DC shift and Normalizing.

Input : input.txt containing the utterances
Output: 1 or 6 for each of utterances

*/

// Config
#define DEFAULT_DURN	"5"
#define N_AMP 			5000 //N_AMP should be high, else it may cause data loss due to scale down - integer round offs
#define N_CLIPS_LIMIT 	15 	// to avoid too many file writes
#define SAMPLERATE 		16000.0
#define ZCR_THR			10 	//%
// Parameters From Initial Analysis
#define M_SILENCE_ENERGY 22500
#define WINDOW_SIZE 		 3790
#define DC_SHIFT 		 0
void writeLines(ofstream &OUTPUT_FS,vector<string> &lines){
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

long long int calcEnergy(vector<int> &samples, int &s, int &e){
	long long int energy= 0;
	for (int i = s; i < e; ++i) {
		energy += (long long)samples[i]*(long long)samples[i];
	}
	return energy;
}
float calcZCR(vector<int> &samples, int &s, int &e){
	int prev= 0,zcr= 0;
	for (int i = s; i < e; ++i) {
		if(samples[i] == 0){
			if(prev * samples[i+1]<0)zcr++;
		}
		else{
			prev=samples[i];//update only if nonzero
			if(samples[i-1] * samples[i]<0)zcr++;
		}
	}
	return 100 * zcr/(float)(e-s);
}
void writeSample(string filename,vector<int> &samples, int &s, int &e,bool verbose=false){
	ofstream OUTPUT_FS(filename);
	OUTPUT_FS<<"SAMPLES: "<<(e-s+1)<<endl;
	OUTPUT_FS<<"BITSPERSAMPLE:	16"<<endl;
	OUTPUT_FS<<"CHANNELS:	1"<<endl;
	OUTPUT_FS<<"SAMPLERATE:	"<<SAMPLERATE<<endl;
	OUTPUT_FS<<"NORMALIZED:	FALSE"<<endl;
	// cout<<s<<" , "<<e<<endl;
	for (int i = s; i <= e; ++i) {
		// cout<<samples[i]<<endl;
		OUTPUT_FS<<samples[i]<<endl;
	}
	OUTPUT_FS.close();
	if(verbose)cout<<"Saved file: "<<filename<<endl;
}








#define FIXED_FLOAT(x) std::fixed <<std::setprecision(3)<<(x) 

int main(int argc, char const *argv[]){
	string file_name="input.txt";
	ifstream INPUT_FS;
	INPUT_FS.open(file_name);

	if(1 || argc>1 || !INPUT_FS.good()){ // .good() returns if file not exists
		string duration=argc>1?argv[1]:DEFAULT_DURN;
		cout<<"Calling recording module, duration:  "<<duration<<"s"<<endl;
		cout<<"Press enter key to start recording: \n";getchar();
		system(("Recording_Module.exe "+duration+" temp.wav "+file_name).c_str());
		// calling open on an already open stream fails//clear flags too
		INPUT_FS.close(); INPUT_FS.clear();
		INPUT_FS.open(file_name);		
		remove("temp.wav");		
	}

	//take file input
	vector <string> lines;
	getLines(lines,INPUT_FS);
	//file loaded in memory, can close it now.
	INPUT_FS.close();

	if(!lines.size()){
		// record instead 
		cout<<"ERROR: Empty data! Please check if correct filename is provided and try again.\n";
		return 1;
	}
	// declarations
	int N_SAMPLES=lines.size()-5;
	vector<int> samples;
	int x,maxX= 0;

	//Read line by line, do DC shift as well as load the samples
	for (int i = 0; i < N_SAMPLES; ++i) {
		sscanf(&lines[i+5][0],"%d",&x);
		x-=DC_SHIFT;
		samples.push_back(x);
		maxX=max(maxX,(int)abs(x));
	}
	//Normalize
	float N_fac=N_AMP/(float)maxX;
	for (int i = 0; i < N_SAMPLES; ++i)
		samples[i] *= N_fac;
	

	vector<bool> blurred_samples(N_SAMPLES,false);
	vector< pair<int,int> > clips;
	// Apply M_SILENCE_ENERGY for windowed threshold on square of amplitude
	long long int windowSum= 0;
	for (int i = 0; i < N_SAMPLES; ++i) {
		x=pow(samples[i],2);//abs(samples[i]);
		windowSum += x;
		if(i>WINDOW_SIZE){
			windowSum -= pow(samples[i-WINDOW_SIZE],2);//abs(samples[i-WINDOW_SIZE]);
			if((windowSum/(float)WINDOW_SIZE)>M_SILENCE_ENERGY)
				blurred_samples[i-WINDOW_SIZE/2]=true;// take midpoint of the window
		}
		// else default is false
	}

	// Detect clips using applied threshold
	bool prev=blurred_samples[0];
	int start= 0;
	for (int i = 1; i < N_SAMPLES; ++i) {
		if(prev != blurred_samples[i]){
			if(!start){
				start=i;
			}else{
				//end
				clips.push_back(make_pair(start,i));
				start= 0;
			}
		}
		prev=blurred_samples[i];
	}
	int clip_size,NUM_CLIPS=clips.size();
	float zcr;double sumZCR= 0;
	long long int energy,sumEnergy= 0;
	if(NUM_CLIPS > N_CLIPS_LIMIT){
		cout<<"Error: too many clips found: "<<NUM_CLIPS<<endl;
		cout<<"Either increase limit or resolve the issue\n";
		return 1;
	}
	else if(NUM_CLIPS== 0){
		cout<<"Error: No clips found: either increase recording volume or check if mic is working.\n";
		return 1;	
	}
	std::ostringstream filename;
	cout<<"Record duration: "<< FIXED_FLOAT(N_SAMPLES/SAMPLERATE)<<"s"<<endl;
	
	char NUMS[2]={'1','6'};
	char DETECTED_NUM='x';//will be updated if detected
	int countNums[2],numIdx;
	for(int i= 0;i<2;i++)countNums[i]= 0;

	// For each clip perform the calculations
		for (int i = 0; i < NUM_CLIPS; ++i) {

			clip_size=(clips[i].second-clips[i].first);
			cout<<"Clip "<<i<<": ";
			cout << "["<< FIXED_FLOAT(clips[i].first/SAMPLERATE)<< 
			","<< FIXED_FLOAT(clips[i].second/SAMPLERATE)<<
			"] \tDuration: "<< FIXED_FLOAT(clip_size/SAMPLERATE)<<"s \t";

			zcr=calcZCR(samples,clips[i].first,clips[i].second);
			cout<<"ZCR = "<< FIXED_FLOAT(zcr)<<"% \t";
			energy=calcEnergy(samples,clips[i].first,clips[i].second);
			cout<<"Energy = "<< energy<<" \t";
			sumZCR+=zcr;
			sumEnergy+=energy;
			if(zcr < ZCR_THR)
				cout <<"= ONE";
			else
				cout <<"= SIX";

			numIdx= (zcr < ZCR_THR)?0:1;
			countNums[numIdx]++;
			if(countNums[numIdx] >= NUM_CLIPS/2)
				DETECTED_NUM=NUMS[numIdx];
			cout<<endl;
		}	
		cout<<"avg ZCR: "<<FIXED_FLOAT(sumZCR/(float)NUM_CLIPS)<<" \tavg Energy: "<<FIXED_FLOAT(sumEnergy/(float)NUM_CLIPS)<<endl;
		if(DETECTED_NUM=='x')
			cout<<"No number was detected more than half times. Hence not saving the clips\n";
		else{
			cout<<"Generating filenames(for saving clips) by overall detected number: "<<DETECTED_NUM<<endl;
			for (int i = 0; i < NUM_CLIPS; ++i) {
			// clear the stream;
				filename.str(std::string());
				filename<<"output_clips/150101021_"<<DETECTED_NUM<<"_"<<i<<".txt";			
			//write cropped sample into new file		
				writeSample(filename.str(),samples,clips[i].first,clips[i].second);
			}
			cout<<"All clips saved in 'output_clips/' directory\n";
		}
		cout<<"Press enter key to close.\n";
		getchar();
		return 0;
	}