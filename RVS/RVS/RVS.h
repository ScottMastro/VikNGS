#pragma once
#include "stdafx.h"
#include "MemoryMapped.h"

#include <string>
#include <vector>
#include <algorithm>

//========================================================
// SNP struct
//========================================================

struct genotypeLikelihood {
	bool valid = false;
	double L00;
	double L01;
	double L11;
};

struct SNP {
	std::string chr;
	int loc;

	std::vector<genotypeLikelihood> gl;
	std::vector<double> p;
	std::vector<double> EG;

	double maf = NULL;

	inline bool operator<(SNP& snp2) { return this->loc < snp2.loc; }
	inline double L00(int index) { return this->gl[index].L00; }
	inline double L01(int index) { return this->gl[index].L01; }
	inline double L11(int index) { return this->gl[index].L11; }

};

inline bool locCompare(SNP lhs, SNP rhs) { return lhs < rhs; }

//========================================================
// functions
//========================================================

//HelperVCF.cpp
std::vector<std::string> parseHeader(MemoryMapped &, int &);
std::vector<bool> getIDs(std::string, std::string, int);
SNP initSNP(MemoryMapped &, std::vector<int>, int);
std::vector<SNP> parseAndFilter(std::string, int, double, std::vector<bool> &);

//Statistics.cpp
std::vector<double> calcEM(SNP &);
std::vector<double> calcEG(SNP &);
std::vector<double> logisticRegression(std::vector<bool> &, std::vector<double> &);
double logisticRegressionInterceptOnly(std::vector<bool> &, std::vector<double> &);
double scoreTest(std::vector<bool> &y, std::vector<double> &x);
double chiSquareOneDOF(double);

//Tests.cpp
std::vector<double> RVSasy(std::vector<SNP> &, std::vector<bool> &, bool );


//========================================================
// inline functions
//========================================================

/**
Removes whitespace from a string

@param str String to remove whitespace from.
@return str without whitespace.
*/
inline std::string trim(std::string str)
{
	str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
	return str;
}

/**
Extracts string from a MemoryMapped class

@param charArray MemoryMapped to extract string from.
@param start Index to start extracting string from.
@param end Index to stop extracting string from.
@return Extracted string.
*/
inline std::string getString(MemoryMapped &charArray, int start, int end) {
	std::string ret;
	for (; start < end; start++) { ret += charArray[start]; }
	return ret;
}

/**
Finds a string in a vector of strings.

@param query The string to find.
@param v The vector to search in.
@return Index of query in v or -1 if query is not found in v.
*/
inline int findIndex(std::string query, std::vector<std::string> v) {
	auto index = std::find(v.begin(), v.end(), query);
	return index != v.end() ? -1 : index - v.begin();
}
