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
#include <assert.h>
using namespace std;

#define ll long long int
#ifdef DEBUG_MODE
#define DEBUG_ON true
#else
#define DEBUG_ON false
#endif
#define DEBUGV(x) if(1 or DEBUG_ON){cout << ">> " << #x << " : \t";for(auto _i = (x).begin();_i!=(x).end();_i++)cout<<(*_i)<<"\t";cout<<"\n";}

/*
Applying K means on universe file
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

double getDistortion(vector< vector<double> > &means,vector<int> &labels, vector< vector<double> > &universe) {
	double dist=0;
	assert(labels.size()==universe.size());
	int N=universe.size();
	for (int i = 0; i < N; ++i){
		dist += vectorDist(means[labels[i]],universe[i])/N;
	}
	return dist;
}
vector< vector<double> >  loadCis(ifstream &INPUT_FS){
	vector< vector<double> >allCis;	
	int frame=0;
	double tempd;
	while(true){
		vector<double>temp(P_ORDER);
		
		// if C[0] is there, skip it
		INPUT_FS>>tempd;

		for (int p = 0; p < P_ORDER; ++p){
			// just like cin>>double_t, it will ignore newlines
			// ignore first coeff
			INPUT_FS>>temp[p];	
		}
		if(INPUT_FS.eof())break;
		allCis.push_back(temp);
	}
	// cout<<filename<<endl;
	return allCis;
}
double rangeFloatRand(double &a, double &b){
	double diff = b-a;
	// RAND_MAX defined in cstdlib header
	return a + diff*rand()/RAND_MAX; 
}
int getLabel(vector<double> &item, vector< vector<double> > &means){
	int label; 
	double dist, minDist = INT_MAX;
	for (int i = 0; i < means.size(); ++i) {
		if(USE_TOKHURA)
			dist = tokhuraDist(means[i],item);
		else
			dist = vectorDist(means[i],item);
		if(dist < minDist){
			minDist=dist;
			label=i;
		}
	}
	return label;
}
void updateMeansAdd(vector<double> &item, vector<double> &means, double &clusterSize){
	// item being added to this mean
	for (int i = 0; i < means.size(); ++i) {
		means[i] = (means[i]*(clusterSize-1)+item[i])/clusterSize;
	}
}
void updateMeansRemove(vector<double> &item, vector<double> &means, double &clusterSize){
	// item being removed from this mean
	for (int i = 0; i < means.size(); ++i) {
		means[i] = (means[i]*(clusterSize+1)-item[i])/clusterSize;
	}
}

vector<vector<double> > k_means(const int K,int iterations, vector< vector<double> >	&universe, double EPSILON=0.30){
	int N = universe.size(), cols=universe[0].size();
	assert(cols==P_ORDER);
	vector< vector<double> > means(K, vector<double>(cols));
	vector<double> clusterSizes(K,0), colmin(cols,(double)INT_MAX), colmax(cols,(double)INT_MIN);
	// calc colwise min maxs
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < cols; ++j){
			colmin[j] = min(colmin[j], universe[i][j]);
			colmax[j] = max(colmax[j], universe[i][j]);
		}
	}
	// initialize means-
	/* For each column, take values randomly between the corresponding minimum and maximum */
	for (int i = 0; i < K; ++i){
		for (int j = 0; j < cols; ++j){
			means[i][j] = rangeFloatRand(colmin[j],colmax[j]);
		}
		// DEBUGV(means[i])
	}
	vector<int> labels(N,-1);
	int newlabel;
	// do iterations
	while(iterations--){
		bool converged=true;
		// make a copy of current means/centroids
		vector< vector<double> > nextMeans = means; 
		for (int i = 0; i < N; ++i)
		{
			// get distance from original means of current iteration
			newlabel = getLabel(universe[i],means);
			if(newlabel!=labels[i]){
				// some label changed => not yet converged
				converged=false;
				if(labels[i]!=-1){
					clusterSizes[labels[i]]--;
					updateMeansRemove(universe[i],nextMeans[labels[i]],clusterSizes[labels[i]]);
				}
				clusterSizes[newlabel]++;
				updateMeansAdd(universe[i],nextMeans[newlabel],clusterSizes[newlabel]);
				labels[i]=newlabel;
			}
		}
		means = nextMeans;
		if(getDistortion(means,labels,universe) < EPSILON)
			break;
		if(converged)
			break;
	}
	for (int i = 0; i < K; ++i){
		DEBUGV(means[i])
	}
	cout << " Labels in order: ";
	DEBUGV(labels)
	cout<< "Remaining iterations: "<<iterations<<endl;
	cout <<"Final distortion: "<< getDistortion(means,labels,universe)<<endl;

	// classify final.
	return means;
}

#define FIXED_FLOAT(x) std::fixed <<std::setprecision(6)<<(x) 
int main(int argc, char const *argv[]) {
	ifstream INPUT_FS;
	INPUT_FS.open("coeffs/universe.txt");
	if(!INPUT_FS.good()){
		cout<<"ERROR: universe not found! xD";
		return 0;
	}
	vector< vector<double> > universe;
	universe = loadCis(INPUT_FS);	
	INPUT_FS.close();
	cout<<"U"<<universe.size()<<endl;
	assert(universe.size());
	// Find the K centroids here-
	vector<vector<double> > codebook = k_means(CODEBOOK_SIZE, KMEANS_ITERATIONS, universe);
	return 0;
}