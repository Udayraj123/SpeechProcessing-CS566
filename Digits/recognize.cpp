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
#define FIXED_FLOAT(x) (x)
// #define FIXED_FLOAT(x) std::fixed <<std::setprecision(6)<<(x) 
set <string> DEBUG_FNs;
#define DEBUG(x) if(DEBUG_ON && DEBUG_FNs.count(__FUNCTION__)){cout << ">> " << #x << " : \t";cout<<FIXED_FLOAT(x)<<endl;}
#define DEBUGV(x) if(DEBUG_ON && DEBUG_FNs.count(__FUNCTION__)){cout << ">> " << #x << " : \t";for(auto _i = (x).begin();_i!=(x).end();_i++)cout<<FIXED_FLOAT(*_i)<<"\t";cout<<"\n";}
#define VAR_NAME(x) #x
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
int NUM_STATES,NUM_SYMBOLS,NUM_OBS;
/*
	Program to generate Parameters used in Hidden Markov Model
	Alpha matrix
	Delta and Psi matrices
*/
// Note: Everything is 0-indexed

string trim(string &aString){
	// https://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
	auto start = aString.find_first_not_of(' ');
	auto end = aString.find_last_not_of(' ');
	return aString.substr(start, (end - start) + 1);
}

void getLines(vector<string>& linesVec_, ifstream& inFile_)
{
	string line;
	while (getline(inFile_, line))
	{
		linesVec_.pb(line);
	}
}
vd  loadSeq(string filename){
	ifstream INPUT_FS;
	INPUT_FS.open(filename);
	
	vd sequence;
	double t;
	while(!INPUT_FS.eof()){
		INPUT_FS>>t;
		sequence.pb(t);
	}
	INPUT_FS.close();
	// cout<<filename<<endl;
	return sequence;
}
void saveMatrix(string filename, vvd &A){
	ofstream OUTPUT_FS(filename);
	for(auto l : A){
		for(auto a : l){
			OUTPUT_FS << a << " ";
		}
		OUTPUT_FS << endl;
	}
	OUTPUT_FS.close();
}
vvd loadMatrix(string filename){
	vector<string> lines;
	ifstream INPUT_FS;
	INPUT_FS.open(filename);
	if(!INPUT_FS.good()){
		cout<<"loadMatrix: File does not exist: "<<filename<<endl;
		exit(1);
	}
	getLines(lines,INPUT_FS);
	INPUT_FS.close();
	
	vvd  matrix;
	double t;
	for(auto &l : lines){
		stringstream ls(l);
		vd temp;
		while(ls>>t)
			temp.pb(t);
		matrix.pb(temp);
	}
	return matrix;
}

void getDeltaPsi(vvd &A,vvd &B,vd &PI,vd &OBS,vvd &Delta,vvi &Psi){
	
	vd row;
	vi psi_row;
	// Initialization
	FOR(i,NUM_STATES){
		row.pb(PI[i]*B[i][OBS[0]]);
		psi_row.pb(0);
	}
	// PRINTV(row)
	Delta[0] = row;
	Psi[0] = psi_row;
	// Induction Step O((T-1)*N^2) 
	FOR(t,NUM_OBS-1){
		vd newrow;
		vi newpsi_row;
		FOR(j,NUM_STATES){
			double maxDel=0;
			int maxI=0;
			FOR(i,NUM_STATES){
				if(maxDel < Delta[t][i]*A[i][j]){
					maxDel = Delta[t][i]*A[i][j];
					maxI = i;
				}
			}
			newrow.pb(maxDel * B[j][OBS[t+1]]);
			newpsi_row.pb(maxI);
		}
		// inserted at t+1
		Delta[t+1] = newrow;
		Psi[t+1] = newpsi_row;
	}
	// FOR(t,NUM_OBS){
		// DEBUGV(Delta[t])
	// }
	// FOR(t,NUM_OBS){
		// DEBUGV(Psi[t])
	// }
}
int main(int argc, char const *argv[]) {
	/** Set visibility of debug messages **/
	DEBUG_FNs.insert("main");
	DEBUG_FNs.insert("saveMatrix");
	// DEBUG_FNs.insert("getDeltaPsi");

	// Label the framewise Cis for each file.
	char DIGITS[] = { '0','1','2','3','4','5','6','7','8','9' };
	int totalMatchCount=0;
	FOR(_UTTERNO,10){
		int matchCount=0;
		for (int _digitIdx = 0; _digitIdx < sizeof(DIGITS)/sizeof(DIGITS[0]); ++_digitIdx)
		{
			char FOUND_DIGIT = 'X';
			int FOUND_UTTER;
			double maxPStar = 0;
			char TEST_DIGIT = DIGITS[_digitIdx];

			std::ostringstream filename2;	
			filename2 << CURR_DIR<< "obs/"<<TEST_DIGIT<<"_"<<_UTTERNO<<"_obs.txt";
			vd OBS = loadSeq(filename2.str());	
			NUM_OBS = min((int)OBS.size(),OBSERVATIONS_LIM);
			// delete excess sequence
			OBS.erase(OBS.begin()+NUM_OBS,OBS.end());
			bool ONE_INDEXED=1;
			for(auto o : OBS) {
				if(!o)
					ONE_INDEXED=0;
			}
			NUM_SYMBOLS=0;	
			for(auto &o : OBS) {
				// make it 0-indexed (just like everything else)
				if(ONE_INDEXED)
					o--;
				NUM_SYMBOLS=max((int)o+1,NUM_SYMBOLS);
			}
			// DEBUG(ONE_INDEXED)

			// PRINT(NUM_OBS)
			// cout<<"Initial Sequence: ";
			// PRINTV(OBS)
			FOR(UTTERNO,10){
				for (int digitIdx = 0; digitIdx < sizeof(DIGITS)/sizeof(DIGITS[0]); ++digitIdx)
				{
					char CHECK_DIGIT = DIGITS[digitIdx];
					std::ostringstream filename;
					filename << CURR_DIR<< "models/"<<CHECK_DIGIT<<"_"<<UTTERNO;

			/** Take Input Data **/
					vd PI = loadSeq("PI.vector");
					vvd A = loadMatrix(filename.str() + "_A.matrix");
					vvd B = loadMatrix(filename.str() + "_B.matrix");
					NUM_STATES=B.size();						
					if(B[0].size() < NUM_SYMBOLS){
						cout<<"B matrix doesn't have enough columns! \n";
						exit(1);
					}				

					double PStar=0; int QStar=0;
					vvd Delta(NUM_OBS,vd(NUM_STATES)); vvi Psi(NUM_OBS,vi(NUM_STATES));
					// Viterbi Algorithm -
					getDeltaPsi(A,B,PI,OBS,Delta,Psi);
					PStar=QStar=0;
					FOR(i,NUM_STATES){
						if(PStar < Delta[NUM_OBS-1][i]){
							PStar = Delta[NUM_OBS-1][i];
							QStar = i;
						}
					}
					// PRINT(NUM_STATES)
					// PRINT(NUM_SYMBOLS)
					// PRINT(PStar)
					// cout<<endl;
					if(maxPStar < PStar){
						maxPStar = PStar;
						FOUND_DIGIT = CHECK_DIGIT;
						FOUND_UTTER = UTTERNO;
					}
				}
			}
			if(TEST_DIGIT==FOUND_DIGIT)matchCount++;
			cout<<(TEST_DIGIT==FOUND_DIGIT?"MATCHED":"NOPE")<<"\t TEST: "<<TEST_DIGIT<<"\t DETECTED: "<<FOUND_DIGIT<<"("<<FOUND_UTTER<<")"<<"\t Max PStar: "<<maxPStar<<endl;

		}
		PRINT(matchCount);
		totalMatchCount += matchCount;
	}
	PRINT(totalMatchCount);
	
	return 0;
}