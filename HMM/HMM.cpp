#include "config.h"
#include <iostream>
#include <string>
#include <vector>
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
#define DEBUG(x) if(DEBUG_ON){cout << ">> " << #x << " : \t";cout<<(x)<<endl;}
#define DEBUGV(x) if(DEBUG_ON){cout << ">> " << #x << " : \t";for(auto i = (x).begin();i!=(x).end();i++)cout<<(*i)<<"\t";cout<<"\n";}
#define PRINT(x) {cout << #x << " : \t";cout<<(x)<<endl;}
#define PRINTV(x){cout << #x << " : \t";for(auto i = (x).begin();i!=(x).end();i++)cout<<(*i)<<" ";cout<<"\n";}
#define FOR(i,n)for (ll i = 0; i < (n); ++i)
#define FORG(i,a,b)for (ll i = (a); i <= (b); ++i)
#define FORD(i,a,b)for (ll i = (a); i >= (b); --i)

#define pb push_back
#define vi vector <int> 
#define vvi vector <vector<int> >
#define vd vector <double> 
#define vvd vector <vector<double> >
int NUM_STATES,SIZE_CODEBOOK,NUM_OBS;

/*
	Program to generate Parameters used in Hidden Markov Model
	Alpha matrix
	Delta and Psi matrices
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
vvd  loadMatrix(string filename){
	vector<string> lines;
	ifstream INPUT_FS;
	INPUT_FS.open(filename);
	getLines(lines,INPUT_FS);
	INPUT_FS.close();
	
	vvd  matrix;
	double t;
	for(auto l : lines){
		stringstream ls(l);
		vd temp;
		while(!ls.eof()){
			ls>>t;
			temp.pb(t);
		}
		matrix.pb(temp);
	}
	// cout<<filename<<endl;
	return matrix;
}


vvd getAlphas(vvd &A,vvd &B,vd &PI,vd &OBS){
	vvd Alpha;

	vd row;
	// Initialization
	FOR(i,NUM_STATES){
		row.pb(PI[i]*B[i][OBS[0]]);
	}
	Alpha.pb(row);
	// Induction Step  O((T-1)*N^2) 
	FOR(t,NUM_OBS-1){
		vd newrow;
		FOR(j,NUM_STATES){
			double sum=0;
			FOR(i,NUM_STATES){
				sum += Alpha[t][i]*A[i][j];
			}
			newrow.pb(sum * B[j][OBS[t+1]]);
		}
		// inserted at t+1
		Alpha.pb(newrow);
	}
	return Alpha;
}
void getDeltaPsi(vvd &A,vvd &B,vd &PI,vd &OBS,vvd &Delta,vvi &Psi){
	
	vd row;
	vi psi_row;
	// Initialization
	FOR(i,NUM_STATES){
		row.pb(PI[i]*B[i][OBS[0]]);
		psi_row.pb(0);
	}
	PRINTV(row)
	Delta.pb(row);
	Psi.pb(psi_row);
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
		PRINTV(newrow)
		// inserted at t+1
		Delta.pb(newrow);
		Psi.pb(newpsi_row);
	}
}

int main(int argc, char const *argv[]) {
	vvd A = loadMatrix("A.matrix");
	vvd B = loadMatrix("B.matrix");
	vd PI = loadSeq("PI.vector");
	vd OBS = loadSeq("O1.txt");	

	NUM_STATES=B.size();
	SIZE_CODEBOOK=B[0].size();
	NUM_OBS = min((int)OBS.size(),OBSERVATIONS_LIM);

	cout<<"Initial Observation Sequence: ";
	FOR(i,NUM_OBS)
		cout<< OBS[i]<<" ";
	cout<<endl;
	// make it 0-indexed
	for(auto &o : OBS)
		o--;

	// Forward Propogation
	vvd Alpha = getAlphas(A,B,PI,OBS);
	// Termination Step	
	double P_Obs_Given_Seq = 0;
	FOR(i,NUM_STATES){
		P_Obs_Given_Seq += Alpha[NUM_OBS-1][i];
	}
	PRINT(P_Obs_Given_Seq);
	
	// Viterbi Algorithm -
	vvd Delta; vvi Psi;
	getDeltaPsi(A,B,PI,OBS,Delta,Psi);
	double PStar=0; int QStar=0;
	FOR(i,NUM_STATES){
		if(PStar < Delta[NUM_OBS-1][i]){
			PStar = Delta[NUM_OBS-1][i];
			QStar = i;
		}
	}
	PRINT(PStar)
	PRINT(QStar)
	vi stateSeq;
	stateSeq.pb(QStar);
	// Backtracking
	FORD(t,NUM_OBS-2,0){
		QStar = Psi[t+1][QStar];
		stateSeq.pb(QStar);
	}
	reverse(stateSeq.begin(),stateSeq.end());
	PRINTV(stateSeq);
	return 0;
}