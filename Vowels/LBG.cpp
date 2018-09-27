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
Applying LBG means on universe file
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
double euclideanDist(vector <double> cis,vector <double> cis2){
	double d = 0 ;
	for (int i = 1; i < cis.size(); ++i){
		d += (cis[i]-cis2[i])*(cis[i]-cis2[i]);
	}
	return d;
}

double getDistortion(vector< vector<double> > &means,vector<int> &labels, vector< vector<double> > &universe) {
	double dist=0;
	assert(labels.size()==universe.size());
	int N=universe.size();
	for (int i = 0; i < N; ++i){
		dist += euclideanDist(means[labels[i]],universe[i])/N;
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
		// if(USE_TOKHURA)
			// dist = tokhuraDist(means[i],item);
		// else
		dist = euclideanDist(means[i],item);
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

vector<vector<double> > LBG(const int K,int iterations, vector< vector<double> >	&universe, double DISTORTION_EPSILON=0.30){
	int N = universe.size(), cols=universe[0].size();
	assert(cols==P_ORDER);
	vector<double> mean_vec(cols,0);
	vector< vector<double> > means;
	vector<double> col_averages(cols,0),split_epsilon(cols);
	// calc colwise min maxs
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < cols; ++j){
			// use as a sum
			col_averages[j] += universe[i][j];
		}
	}
	double sum_dev;
	for (int j = 0; j < cols; ++j){
		col_averages[j]/=N;
		sum_dev=0;
		for (int i = 0; i < N; ++i){
			sum_dev += (universe[i][j]-col_averages[j])*(universe[i][j]-col_averages[j]);
		}		
		// std deviation / 1000
		split_epsilon[j] = sqrt(sum_dev/N)/1000;
		// split_epsilon[j] = sqrt(sum_dev/N)/100;
	}
	// DEBUGV(split_epsilon)

	// start with 1 vector = colwise mean of universe
	means.push_back(col_averages);

	vector<int> labels(N,-1);			
	// do iterations
	while(means.size()<K){
		// make a copy of current means/centroids
		vector< vector<double> > nextMeans = means;
		// Do the splitting for every centroid 
		for (int i = 0; i < means.size(); ++i){
			for (int j = 0; j < cols; ++j){
				nextMeans[i][j] -= split_epsilon[j];
				means[i][j] += split_epsilon[j];
			}
			// DEBUGV(means[i])
			// DEBUGV(nextMeans[i])
		}
		// double the size
		means.insert(means.end(),nextMeans.begin(),nextMeans.end());
		// use as a sum
		vector<vector<double> > nextColSums(means.size(),std::vector<double>(cols,0));
		for (int i = 0; i < N; ++i)
		{
			// get distance from original means of current iteration
			labels[i] = getLabel(universe[i],means);
			for (int j = 0; j < cols; ++j)
			{
				nextColSums[labels[i]][j] += universe[i][j];
			}
		}

		// update to new col wise averages for each centroid
		for (int i = 0; i < nextColSums.size(); ++i){
			for (int j = 0; j < cols; ++j){
				means[i][j]= nextColSums[i][j]/N;
			}
		}
		if(getDistortion(means,labels,universe) < DISTORTION_EPSILON)
			break;
	}

	for (int i = 0; i < CODEBOOK_SIZE; ++i) {
		cout<<"i=" <<i <<" ";DEBUGV(means[i]);
	}
	cout << " Labels in order: ";
	DEBUGV(labels);
	cout<< "Final codebook size: "<<means.size()<<endl;
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
	vector<vector<double> > codebook = LBG(CODEBOOK_SIZE, KMEANS_ITERATIONS, universe);	
	return 0;
}