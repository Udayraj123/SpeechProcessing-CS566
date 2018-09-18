// #include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

constexpr auto DC = 0;		//DC Offset (Obtained to be approx. 0 for my setup);
constexpr auto ZT = 10;		//ZCR Threshhold percentage(Obtained Emperically using debug statements);
using namespace std;

std::vector <double> hamming(std::vector<double> s, int wsize = -1)
{
	if (wsize < 0)
		wsize = s.size();
	for (int i = 0; i < wsize; i++)
	{
		s[i] = s[i]*(0.56 - 0.46*cos((2 * M_PI * i) / wsize));
	}
	return s;
}

std::vector<double> calc_R(std::vector<double> s, int p)
{
	std::vector<double> R(p+1,0);
	for (int i = 0; i <= p; i++)
	{
		for (int m = 0; i + m < s.size(); m++)
			R[i] += s[m]*s[i + m];
	}
	return R;
}

std::vector<double> calc_A(std::vector<double> &R)
{
	int p = R.size();
	double k;
	std::vector<double> A(p, 0), E(p, 0), A_temp;
	E[0] = R[0];

	for (int i = 1; i < p; i++)
	{
		A_temp = A;
		k=R[i];
		for (int j = 1; j < i; j++)
			k -= A_temp[j] * R[i - j];
		k /= E[i - 1];
		A[i] = k;
		for (int j = 1; j < i; j++)
			A[j] = A_temp[j] - k * A_temp[i - j];
		E[i] = (1 - k * k)*E[i - 1];
	}
	return A;
}

std::vector<double> calc_C(std::vector<double> &A, double sigma_sq)
{
	int p = A.size();
	std::vector<double> C(A);
	C[0] = log2(sigma_sq);
	for (int i = 1; i < p; i++)
	{
		for (int m = 1; m < i; m++)
			C[i] += m * C[m] * A[i - m] / i;
	}
	return C;
}

double dist(std::vector<double> ref, std::vector<double> target, std::vector<double> weight)
{
	double res = 0, diff;
	for (int i=1; i<ref.size(); i++)
	{
		diff = ref[i] - target[i];
		res += weight[i] * diff*diff;
	}
	return res;
}

std::vector <double> readFile(std::string inputfile)
{
	std::ifstream fin;
	std::vector <double> v(0);
	//DEBUG STATEMENT
	//inputfile = "D:\\Programs\\Cooledit\\150101012_" + std::to_string(num) + "_" + std::to_string(file_num) + ".txt";

	fin.open(inputfile);

	if (fin.fail())
	{
		std::cout << "Unable to open \"" << inputfile << "\"!\n\n";
		std::cin.ignore(1000, '\n');
		return v;
	}

	std::string tempstr;

	while (true)				//SKIP INITIAL LINES
	{
		std::getline(fin, tempstr);
		if (tempstr[0] == '-' || tempstr[0] < 58 && tempstr[0]>47)
			break;
	}

	double inp ;//= stod(tempstr);
	sscanf(&tempstr[0], "%lf", &inp);
	int max = 1;
	while (!fin.eof())			//READ SAMPLES
	{
		inp -= DC;				//DC shift
		v.push_back(inp);
		if (abs(inp) > max)
			max = inp;
		fin >> inp;
	}
	double nrm = (double)5000 / max;	//Normalisation factor
	for (int i = 0; i < v.size(); i++)
	{
		v[i] *= nrm;
	}
	fin.close();
	return v;
}

std::pair<int, int> VAD(std::vector <double> v)
{
	double e = 0;
	int i;
	for (i = 0; i < 400; i++)	//CALCULATE ENERGY RANGE FOR SILENCE(NOISE)
	{
		e += v[i] * v[i];
	}
	double et = 625 * e;			//SET ENERGY THRESHHOLD

	for (; i < v.size(); i++)	//FIND START POINT
	{
		e += v[i] * v[i] - v[i - 400] * v[i - 400];
		if (e > et)
			break;
	}
	int start = i - 400;
	int end = v.size();
	for (; i < v.size(); i++)	//FIND END POINT
	{
		e += v[i] * v[i] - v[i - 400] * v[i - 400];
		if (e < et)
		{
			end = i - 400;
			break;
		}
	}
	return std::make_pair(start, end);
}

int main()
{
	std::ofstream fout;

	//DEBUG STATEMENT
	//int num = 1, file_num=1;

	std::cout << "NOTE: Enter full file path if file is not in present present working directory!\n";
	char sel = 'y';
	do
	{
		//DEBUG STATEMENT
		//inputfile = "D:\\Programs\\Cooledit\\150101012_" + std::to_string(num) + "_" + std::to_string(file_num) + ".txt
		std::string inputfile;
		std::cout << "Input file (No Spaces) : ";
		std::cin >> inputfile;
		std::vector <double> v = readFile(inputfile);
		std::vector <double> w{ 1,3,7,13,19,22,25,33,42,50,56,61};
		if (v.size() < 1)
			continue;

		std::string outputfile="Ci.txt";
		fout.open(outputfile);

		int p = 12, window_size = 320, stride = 80, no_of_itern=5;
		std::vector <double> R, A, C;

		int start=0, end=v.size();
		std::pair<int, int> tempp = VAD(v);
		start=tempp.first;
		end=tempp.second;
		start += (end - start - window_size - stride*(no_of_itern-1) - 1) / 2;
		start = v.size()/4;
		for (int itern=0;itern<no_of_itern;itern++)
		{
			std::vector <double> samples(v.begin() + start + stride * itern, v.begin() + start + stride * itern + window_size);
			
			//R = calc_R(samples, p);
			R = calc_R(hamming(samples), p);
			std::cout << "\nIteration " << itern << ":\nR=[";
			for (int i = 0; i <= p; i++)
			{
				std::cout << R[i]/R[0]<<' ';
			}

			A = calc_A(R);
			std::cout << "]\nA=[";
			for (int i = 1; i <= p; i++)
			{
				std::cout << A[i] << ' ';
				//fout << A[i] << ' ';
			}
			//fout << '\n';

			C = calc_C(A,R[0]);
			std::cout << "]\nC=[";
			for (int i = 0; i <= p; i++)
			{
				std::cout << C[i] << ' ';
				fout << C[i] << ' ';
			}
			std::cout << "]\n";
			fout << '\n';
			/*for (;;)
			{
				Tdist[i]=dist(ref[i], v, w)
			}*/
		}
		std::cout << "\nRun Again?(y/N) : ";
		std::cin >> sel;
		std::cin.ignore();
		fout.close();
	} while (sel == 'y' || sel == 'Y');
	return 0;
}

/*
//Calculate ZCR
bool lastsign = (v[start] > 0);
int z = 0;
for (i = start + 1; i < end; i++)
{
	if (v[i] == 0)
		continue;
	if ((v[i] > 0) ^ lastsign)
		lastsign = !lastsign, z++;
}

//DEBUG STATEMENT
// std::cout << inputfile << '\t' << v.size() << '\t' << start << '\t' << end << '\t' << z * 100 / (end - start) << '\n';

if (100 * z < ZT * (end - start))	//ZCR COMPARISION
	std::cout << "ONE\n";
else
	std::cout << "SIX\n";
fin.close();

//DEBUG STATEMENTS
//	file_num++;
//	if (file_num > 10)
//		num = 6, file_num = 1;
//} while (!(num==6 && file_num==10));
*/
