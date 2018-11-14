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


vvd getAlpha(vvd &A,vvd &B,vd &PI,vd &OBS){
	vvd Alpha(NUM_OBS,vd(NUM_STATES));
	vd row;
	// Initialization
	FOR(i,NUM_STATES){
		row.pb(PI[i]*B[i][OBS[0]]);
	}
	Alpha[0] = row;
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
		Alpha[t+1] = newrow;
	}
	FOR(t,NUM_OBS){
		DEBUGV(Alpha[t])
	}

	return Alpha;
}
vvd getBeta(vvd &A,vvd &B,vd &OBS){
	vvd Beta(NUM_OBS,vd(NUM_STATES));

	// Initialization
	vd row(NUM_STATES,1);
	Beta[NUM_OBS-1]=row;

	// Induction Step  O((T-1)*N^2) 
	FORD(t,NUM_OBS-2,0){
		// Note: Everything is 0-indexed
		vd newrow;
		FOR(i,NUM_STATES){
			double sum=0;
			FOR(j,NUM_STATES){
				sum += A[i][j]*Beta[t+1][j]*B[j][OBS[t+1]];
			}
			newrow.pb(sum);
		}
		Beta[t]=newrow;
	}
	FOR(t,NUM_OBS){
		DEBUGV(Beta[t])
	}
	return Beta;
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
vvvd getXImatrix3D(vvd &A,vvd &B,vd &OBS,vvd &Alpha,vvd &Beta){
	vvvd XImatrix3D(NUM_OBS,vvd(NUM_STATES,vd(NUM_STATES)));
	double sum;
	FOR(t,NUM_OBS-1){
		sum=0;
		FOR(i, NUM_STATES){
			FOR(j, NUM_STATES){
				sum += (XImatrix3D[t][i][j] = Alpha[t][i]*A[i][j]*B[j][OBS[t+1]]*Beta[t+1][j]);
			}
		}
		DEBUG(t)
		FOR(i,NUM_STATES){
			FOR(j,NUM_STATES){
				// XImatrix3D[t][i][j] /= sum;
				XImatrix3D[t][i][j] = sum ? (XImatrix3D[t][i][j]/sum): 0;
			}
			DEBUGV(XImatrix3D[t][i])
		}
	}
	return XImatrix3D;
}

vvd getGamma(vvvd &XImatrix3D){
	vvd Gamma(NUM_OBS,vd(NUM_STATES));
	FOR(t,NUM_OBS-1){
		FOR(i,NUM_STATES){
			double sum = 0;
			FOR(j,NUM_STATES){
				sum+=XImatrix3D[t][i][j];
			}
			Gamma[t][i]=sum;
		}
		if(t<5){
			DEBUG(t)
			DEBUGV(Gamma[t])
		}
	}
	return Gamma;
}
void updateABPI(vvd &A,vvd &B,vd &PI,vd &OBS,vvvd &XImatrix3D,vvd &Gamma){
	double gSum,sum;

	// Update AB
	FOR(i,NUM_STATES){
		gSum=0;
		vd symSum(NUM_SYMBOLS,0);
		FOR(t,NUM_OBS-1){
			// For A
			gSum+=Gamma[t][i];
			// For B
			symSum[OBS[t]] += Gamma[t][i];
		}
		// Update A
		FOR(j,NUM_STATES){
			sum=0;
			FOR(t,NUM_OBS-1){
				sum+=XImatrix3D[t][i][j];
			}
			A[i][j]=sum/gSum;
		}
		// Update B
		FOR(k,NUM_SYMBOLS){
			B[i][k]=symSum[k]/gSum;
		}
	}
	// Update PI
	FOR(i,NUM_STATES){
		PI[i]=Gamma[0][i];
	}
}

template<class T>
void operator /= (std::vector< vector<T> >& vec, const T& N) {
	for(auto &l : vec){
		for(auto &v : l){
			v /= N;
		}
	}
};
template<class T>
void operator += (std::vector< vector<T> >& v1, std::vector< vector<T> >& v2) {
	int i=0;
	for(auto &l : v1){
		int j=0;
		for(auto &v : l){
			v += v2[i][j];
			j++;
		}
		i++;
	}
};

int main(int argc, char const *argv[]) {
	/** Set visibility of debug messages **/
	DEBUG_FNs.insert("main");
	DEBUG_FNs.insert("saveMatrix");
	// DEBUG_FNs.insert("getAlpha");
	// DEBUG_FNs.insert("getBeta");
	// DEBUG_FNs.insert("getDeltaPsi");
	// DEBUG_FNs.insert("getXImatrix3D");
	// DEBUG_FNs.insert("getGamma");
	// DEBUG_FNs.insert("updateABPI");

	// Label the framewise Cis for each file.
	char DIGITS[] = { '0','1','2','3','4','5','6','7','8','9' };
	std::ostringstream filename;
	vvd A_init = loadMatrix("A.matrix");
	vvd B_init = loadMatrix("B.matrix");
	vd PI_init = loadSeq("PI.vector");

	FOR(digitIdx,sizeof(DIGITS)/sizeof(DIGITS[0])){
		char TRAIN_DIGIT = DIGITS[digitIdx];
		FOR(EPOCHNO, NUM_EPOCHS){
			vvd A_sum(A_init.size(),vd(A_init[0].size(),0));
			vvd B_sum(B_init.size(),vd(B_init[0].size(),0));

			filename.str(std::string());
			filename << CURR_DIR<< "models/"<<TRAIN_DIGIT<<"_e"<<EPOCHNO;
			
			vvd A_curr = EPOCHNO ? loadMatrix(filename.str()+"_A.matrix") : A_init;
			vvd B_curr = EPOCHNO ? loadMatrix(filename.str()+"_B.matrix") : B_init;

			FOR(UTTERNO,10){
				filename.str(std::string());
				filename << CURR_DIR<< "obs/"<<TRAIN_DIGIT<<"_"<<UTTERNO<<"_obs.txt";
			/** Take Input Data **/
				vd OBS = loadSeq(filename.str());	
				
				vvd A = A_curr;
				vvd B = B_curr;
				vd PI =	PI_init;
			/** Preprocess Input Data **/
				NUM_OBS = min((int)OBS.size(),OBSERVATIONS_LIM);
				// delete excess sequence
				OBS.erase(OBS.begin()+NUM_OBS,OBS.end());
				NUM_STATES=B.size();	

				bool ONE_INDEXED=1;
				for(auto o : OBS) {
					if(!o)
						ONE_INDEXED=0;
				}
				// DEBUG(ONE_INDEXED)
				NUM_SYMBOLS=0;	
				for(auto &o : OBS) {
					// make it 0-indexed (just like everything else)
					if(ONE_INDEXED)
						o--;
					NUM_SYMBOLS=max((int)o+1,NUM_SYMBOLS);
				}
				if(B[0].size() < NUM_SYMBOLS){
					cout<<"B matrix doesn't have enough columns! \n";
					exit(1);
				}

				// PRINT(NUM_STATES)
				// PRINT(NUM_OBS)
				// PRINT(NUM_SYMBOLS)

				// cout<<"Initial Sequence: ";
				// PRINTV(OBS)
				// cout<<endl;

			/** Solution to Problem 1 & 2 **/

				// Forward Procedure
				vvd Alpha = getAlpha(A,B,PI,OBS);
				// Termination Step	
				double P_Obs_Given_Seq = 0;
				FOR(i,NUM_STATES){
					P_Obs_Given_Seq += Alpha[NUM_OBS-1][i];
				}

				// PRINT(P_Obs_Given_Seq);	
				// Backward Procedure
				vvd Beta = getBeta(A,B,OBS);

			/** Solution to Problem 3 **/

				double prevPStar, PStar=0; int QStar=0;
				vvd Delta(NUM_OBS,vd(NUM_STATES)); vvi Psi(NUM_OBS,vi(NUM_STATES));
				vvvd XImatrix3D;
				vvd Gamma;

				FOR(iter,NUM_ITERATIONS){		
					PRINT(iter)
					if(iter>0){
						// dont update for first iter
						XImatrix3D = getXImatrix3D(A,B,OBS,Alpha,Beta);
						Gamma = getGamma(XImatrix3D);
						updateABPI(A,B,PI,OBS,XImatrix3D,Gamma);
					}

					// Viterbi Algorithm -
					getDeltaPsi(A,B,PI,OBS,Delta,Psi);
					PStar=QStar=0;
					FOR(i,NUM_STATES){
						if(PStar < Delta[NUM_OBS-1][i]){
							PStar = Delta[NUM_OBS-1][i];
							QStar = i;
						}
					}
					PRINT(PStar)
					if(iter>0)
						PRINT(PStar - prevPStar)
					// PRINT(QStar)
					cout<<endl;

					prevPStar = PStar;
				}
				if(EPOCHNO == NUM_EPOCHS-1){
					filename.str(std::string());
					filename << CURR_DIR<< "models/"<<TRAIN_DIGIT<<"_"<<UTTERNO;
					// save the matrices as part of model.
					saveMatrix(filename.str()+"_A.matrix", A);
					saveMatrix(filename.str()+"_B.matrix", B);
				}
				// overloaded operator
				A_sum += A;
				B_sum += B;
				// FOR(i,NUM_STATES){
					// DEBUGV(A[i])
				// }
				// cout<<endl;
				// FOR(i,NUM_STATES){
					// DEBUGV(B[i])
				// }
				// cout<<endl;

				// vi FinalStateSeq(NUM_OBS);
				// Make states back to 1 index: 
				// FinalStateSeq[NUM_OBS-1] = (QStar + !ONE_INDEXED);
				// Backtracking
				// FORD(t,NUM_OBS-2,0){
				// QStar = Psi[t+1][QStar];
				// FinalStateSeq[t] = (QStar+ !ONE_INDEXED);
			// }
			// PRINTV(FinalStateSeq);		
				
			}
			// average of A matrix over the utterances will be the start point for new.
			A_sum/=10.0;
			B_sum/=10.0;

			filename.str(std::string());
			filename << CURR_DIR<< "models/"<<TRAIN_DIGIT<<"_e"<<(EPOCHNO+1);
			saveMatrix(filename.str()+"_A.matrix", A_sum);
			saveMatrix(filename.str()+"_B.matrix", B_sum);
		}
	}
	return 0;
}