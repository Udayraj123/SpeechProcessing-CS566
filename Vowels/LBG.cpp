#include "config.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
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
#define FIXED_FLOAT(x) std::fixed <<std::setprecision(4)<<(x) 
#define DEBUG(x) if(1 or DEBUG_ON){cout << ">> " << #x << " : \t"<<FIXED_FLOAT(x)<<endl;}
#define DEBUGV(x) if(1 or DEBUG_ON){cout << ">> " << #x << " : \t";for(auto _i = (x).begin();_i!=(x).end();_i++)cout<<FIXED_FLOAT(*_i)<<"\t";cout<<"\n";}

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
double tokhuraDist(vector <double> &cis,vector <double> &cis2){
	static double weights[]={1,3,7,13,19,22,25,33,42,50,56,61};
	double d = 0 ;
	for (int i = (cis.size()>P_ORDER); i < cis.size(); ++i){		
		d += (cis[i]-cis2[i])*(cis[i]-cis2[i])*weights[i-1];
	}
	return d;
}
double euclideanDistSq(vector <double> &cis,vector <double> &cis2){
	double d = 0 ;
	for (int i = (cis.size()>P_ORDER); i < cis.size(); ++i){
		d += (cis[i]-cis2[i])*(cis[i]-cis2[i]);
	}
	return d;
}

double getDistortion(vector< vector<double> > &codebook,vector<int> &labels, vector< vector<double> > &universe) {
	assert(labels.size()==universe.size());
	double dist=0;
	int N=universe.size();
	for (int i = 0; i < N; ++i){
		dist += euclideanDistSq(codebook[labels[i]],universe[i])/N;
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
	// cout<<filename<<endl;
	return allCis;
}
double rangeFloatRand(double a, double b){
	double diff = b-a;
	// RAND_MAX defined in cstdlib header
	return a + diff*rand()/RAND_MAX; 
}
int getLabel(vector<double> &item, vector< vector<double> > &codebook){
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

vector<double> calcSplit(vector<double> &mean, vector<int> &labels, int label, vector< vector<double> >	&universe){
	int N = universe.size(), cols=universe[0].size();
	double sum_dev;
	vector<double> split_epsilon(cols);
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

vector<vector<double> > handleEmptyCells(vector<vector<double> > &codebook, vector<double> &is_empty, vector<int> &labels){
	//Deal with Empty cell problem
	DEBUGV(is_empty)
	vector<vector<double> > newMeans;
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
vector<vector<double> > LBG(int K,int iterations, vector< vector<double> >	&universe, double DISTORTION_EPSILON=0.1){
	K = pow(2,ceil(log2(K)));	
	const int N = universe.size(), cols=universe[0].size();
	assert(cols==P_ORDER);
	vector<double> col_averages(cols,0);
	vector< vector<double> > codebook;
	
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
	vector<int> labels(N,0);
	// splitting parameter dynamically calculated			
	vector<double> split_epsilon;
	
	double minDistortion=INT_MAX, currDistortion;
	// do iterations
	while(codebook.size()<K){
		// make a copy of current codebook/centroids
		vector< vector<double> > nextCodebook(codebook.begin(),codebook.end());
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
			// std::vector<bool> is_empty(codebook.size(),1);
			nochange=true;
			for (int i = 0; i < N; ++i)
			{
				newlabel = getLabel(universe[i],codebook);
				nochange = nochange && (labels[i]==newlabel);
				labels[i] = newlabel;
				// is_empty[newlabel] = 0;
			}
			// codebook = handleEmptyCells(codebook, is_empty, labels);
			vector<vector<double> > nextColSums(codebook.size(),std::vector<double>(cols,0));
			vector<double> labelCount(codebook.size(),0);
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
			
			if(DEBUG_ON)
				cout<<"iter "<<iter<<": "<<currDistortion<<" / "<<DISTORTION_EPSILON<<endl;

			if(nochange){
				if(DEBUG_ON)
					cout<<"No further label change\n";
				break;
			}
		}
		while(0<--iter && currDistortion > DISTORTION_EPSILON);
	}
	cout<<"Final centroids: \n";
	for (int i = 0; i < CODEBOOK_SIZE; ++i) {
		cout<<"i=" <<i <<" ";DEBUGV(codebook[i]);
	}
	cout << " Labels in order (multiple labels per vowel if K>5): ";
	DEBUGV(labels);
	cout<< "Final codebook size: "<<codebook.size()<<endl;
	cout <<"Final distortion: "<< getDistortion(codebook,labels,universe)<<endl;
	cout <<"Minimum distortion: "<< minDistortion<<endl;
	// classify final.
	return codebook;
}


int main(int argc, char const *argv[]) {
	ifstream INPUT_FS;
	INPUT_FS.open("coeffs/universe.txt");
	if(!INPUT_FS.good()){
		cout<<"ERROR: universe not found! xD";
		return 0;
	}
	vector< vector<double> > universe = loadCis(INPUT_FS);	
	INPUT_FS.close();
	DEBUG(universe.size())
	assert(universe.size());
	// Find the K centroids here-
	vector<vector<double> > codebook = LBG(CODEBOOK_SIZE, KMEANS_ITERATIONS, universe);	
	return 0;
}