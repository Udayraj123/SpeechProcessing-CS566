#include "config.h"
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set> // debug fn names list
#include <algorithm> //reverse and sort
#include <fstream>
#include <sstream> // for i in filename
#include <iomanip> // for setprecision
# define M_PI           3.14159265358979323846
#include <math.h> //for log2,abs()
#include <stdlib.h> // for system()
#include <assert.h>
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
Applying LBG codebook on universe file
*/
void getLines(vector<string>& linesVec_, ifstream& inFile_)
{
	string line;
	while (getline(inFile_, line))
	{
		linesVec_.push_back(line);
	}
}
double tokhuraDist(vd &cis,vd &cis2){
	static double weights[]={1,3,7,13,19,22,25,33,42,50,56,61};
	double d = 0 ;
	for (int i = (cis.size()>P_ORDER); i < cis.size(); ++i){		
		d += (cis[i]-cis2[i])*(cis[i]-cis2[i])*weights[i-1];
	}
	return d;
}
double euclideanDistSq(vd &cis,vd &cis2){
	double d = 0 ;
	for (int i = (cis.size()>P_ORDER); i < cis.size(); ++i){
		d += (cis[i]-cis2[i])*(cis[i]-cis2[i]);
	}
	return d;
}

double getDistortion(vvd &codebook,vi &labels, vvd &universe) {
	assert(labels.size()==universe.size());
	double dist=0;
	int N=universe.size();
	for (int i = 0; i < N; ++i){
		dist += euclideanDistSq(codebook[labels[i]],universe[i])/N;
	}
	return dist;
}
vvd loadCis(string filename){
	ifstream INPUT_FS;
	INPUT_FS.open(filename);
	if(!INPUT_FS.good()){
		cerr<<"loadCis: File does not exist: "<<filename<<endl;
		exit(1);
	}
	vvd allCis;	
	int frame=0;
	double tempd;
	while(true){
		vd temp(P_ORDER);
		
		// if C[0] is there, skip it
		if(SKIP_FIRST_CI)
			INPUT_FS>>tempd;

		for (int p = 0; p < P_ORDER; ++p){
			// just like cin>>double_t, it will ignore newlines
			// ignore first coeff
			INPUT_FS>>temp[p];	
		}
		if(INPUT_FS.eof())break;
		allCis.push_back(temp);
	}
	INPUT_FS.close();

	// cout<<filename<<endl;
	return allCis;
}
double rangeFloatRand(double a, double b){
	double diff = b-a;
	// RAND_MAX defined in cstdlib header
	return a + diff*rand()/RAND_MAX; 
}
int getLabel(vd &item, vvd &codebook){
	int label=-1; 
	double dist, minDist = INT_MAX;//std::numeric_limits<double>::max();
	for (int i = 0; i < codebook.size(); ++i) {
		if(USE_TOKHURA)
			dist = tokhuraDist(codebook[i],item);
		else
			dist = euclideanDistSq(codebook[i],item);
		if(dist < minDist){
			minDist=dist;
			label=i;
		}
	}
	return label;
}

vd calcSplit(vd &mean, vi &labels, int label, vvd &universe){
	int N = universe.size(), cols=universe[0].size();
	double sum_dev;
	vd split_epsilon(cols);
	for (int j = 0; j < cols; ++j){
		sum_dev=0;
		for (int i = 0; i < N; ++i){
			if(labels[i]==label)
				sum_dev += pow(universe[i][j]-mean[j],2);
		}		
		// std deviation / 1000
		split_epsilon[j] = sqrt(sum_dev/N)/1000.0;
		// split_epsilon[j] = rangeFloatRand(0.001,0.002);
	}
	return split_epsilon;
}

vvd handleEmptyCells(vvd &codebook, vd &is_empty, vi &labels){
	//Deal with Empty cell problem
	DEBUGV(is_empty)
	vvd newMeans;
	map<int,int>labelMap;
	int counter=0;
	for (int i = 0; i < codebook.size(); ++i){
		if(!is_empty[i]){
			newMeans.push_back(codebook[i]);								
			labelMap[i]=counter++;
		}
	}
	// for(auto k : labelMap) cout<< k.first <<" : "<<k.second<<", ";  cout<< endl;
	// Shift label numbers when empty cell is removed
	for (int i = 0; i < labels.size(); ++i) {
		assert(labelMap.find(labels[i])!=labelMap.end());
		labels[i]=labelMap[labels[i]];
	}
	return newMeans;
}
vvd LBG(int K,int iterations, vvd &universe, double DISTORTION_EPSILON=0.1){
	K = pow(2,ceil(log2(K)));	
	const int N = universe.size(), cols=universe[0].size();
	assert(cols==P_ORDER);
	vd col_averages(cols,0);
	vvd codebook;
	
	// calc colwise averages
	for (int j = 0; j < cols; ++j){
		col_averages[j]=0;
		for (int i = 0; i < N; ++i){
			// use as a sum
			col_averages[j] += universe[i][j];
		}
		col_averages[j]/=N;
	}
	// start with 1 vector = colwise mean of universe
	codebook.push_back(col_averages);
	// with their label as 0
	vi labels(N,0);
	// splitting parameter dynamically calculated			
	vd split_epsilon;
	
	double minDistortion=INT_MAX, currDistortion;
	// do iterations
	while(codebook.size()<K){
		// make a copy of current codebook/centroids
		vvd nextCodebook(codebook.begin(),codebook.end());
		// Do the splitting for every centroid 
		for (int i = 0; i < codebook.size(); ++i){
			split_epsilon = calcSplit(codebook[i],labels,i,universe);
			for (int j = 0; j < cols; ++j){
				nextCodebook[i][j] -= split_epsilon[j];
				codebook[i][j] += split_epsilon[j];
			}
			// DEBUGV(split_epsilon)
			// DEBUGV(codebook[i])
			// DEBUGV(nextCodebook[i])
		}
		// double the size
		codebook.insert(codebook.end(),nextCodebook.begin(),nextCodebook.end());
		DEBUG(codebook.size())
		bool nochange;
		int iter = iterations, newlabel;
		do
		{
			//  vector<bool> is_empty(codebook.size(),1);
			nochange=true;
			for (int i = 0; i < N; ++i)
			{
				newlabel = getLabel(universe[i],codebook);
				nochange = nochange && (labels[i]==newlabel);
				labels[i] = newlabel;
				// is_empty[newlabel] = 0;
			}
			// codebook = handleEmptyCells(codebook, is_empty, labels);
			vvd nextColSums(codebook.size(), vd(cols,0));
			vd labelCount(codebook.size(),0);
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < cols; ++j){
					nextColSums[labels[i]][j] += universe[i][j];
				}
				labelCount[labels[i]]++;
			}
		// update to new col wise averages for each centroid
			for (int i = 0; i < codebook.size(); ++i){
				if(labelCount[i]==0){
					cout<<"WARNING: Empty cell problem detected!\n";
					continue;
				}
				for (int j = 0; j < cols; ++j){
					// SUCH A BLUNDER WITH THE N!
					// codebook[i][j]= nextColSums[i][j]/N;
					codebook[i][j]= nextColSums[i][j]/labelCount[i];
				}
			}						
			currDistortion = getDistortion(codebook,labels,universe);
			minDistortion = min(minDistortion,currDistortion);
			
			// if(DEBUG_ON)
			// 	cout<<"iter "<<iter<<": "<<currDistortion<<" / "<<DISTORTION_EPSILON<<endl;

			if(nochange){
				// DEBUG("No further label change\n");
				break;
			}
		}
		while(0<--iter && currDistortion > DISTORTION_EPSILON);
	}
	cout<<"Final centroids: \n";
	for (int i = 0; i < CODEBOOK_SIZE; ++i) {
		cout<<"i=" <<i <<" ";DEBUGV(codebook[i]);
	}
	cout << " Labels in order: ";
	DEBUGV(labels);
	cout<< "Final codebook size: "<<codebook.size()<<endl;
	cout <<"Final distortion: "<< getDistortion(codebook,labels,universe)<<endl;
	cout <<"Minimum distortion: "<< minDistortion<<endl;
	// classify final.
	return codebook;
}

// Ad-hoc polymorphism is exhibited by overloading in C++
// Coercion polymorphism is aka implicit Type casting.
// >> For a Class one can define what to do on implicit type casting by operator overloading, here using the other type say int() operator will be overloaded
// Parametric polymorphism aka compile-time polymorphism is exhibited by templates in C++
template <class Vec> void writeVec(string filename, Vec &samples, int s=0, int e=-1, bool verbose = false) {
if(e==-1)e=samples.size()-1;
ofstream OUTPUT_FS(filename);
for (int i = s; i < e; ++i) {
	OUTPUT_FS << samples[i] << " ";
}
OUTPUT_FS << samples[e]<<endl; 
OUTPUT_FS.close();
if (verbose)cout << "Saved file: " << filename << endl;
}

int main(int argc, char const *argv[]) {
	DEBUG_FNs.insert("main");
	DEBUG_FNs.insert("LBG");
	ifstream INPUT_FS;
	vvd universe = loadCis("coeffs/universe.txt");	
	DEBUG(universe.size())
	// assert(universe.size());
	// // Find the K centroids here-
	vvd codebook = LBG(CODEBOOK_SIZE, KMEANS_ITERATIONS, universe);	

	// Label the framewise Cis for each file.
	char DIGITS[] = { '0','1','2','3','4','5','6','7','8','9' };
	std::ostringstream filename;
	FOR(UTTERNO,10){
		for (int digitIdx = 0; digitIdx < sizeof(DIGITS)/sizeof(DIGITS[0]); ++digitIdx)
		{
			char TRAIN_DIGIT = DIGITS[digitIdx];
			filename.str(std::string());
			filename << CURR_DIR<< "coeffs/"<<TRAIN_DIGIT<<"_"<<UTTERNO<<"_cis.txt";
			ifstream INPUT_FS;			
			vvd digitCis = loadCis(filename.str());
			vi obsSeq;
			for(auto v : digitCis){
				obsSeq.push_back(getLabel(v,codebook));
			}
			// DEBUGV(obsSeq)
			filename.str(std::string());
			filename << CURR_DIR<< "obs/"<<TRAIN_DIGIT<<"_"<<UTTERNO<<"_obs.txt";
			writeVec(filename.str(),obsSeq);
		}
	}
	return 0;
}