#include "stdafx.h"
#include "RVS.h"
#include "MemoryMapped.h"

#include <iostream>
#include <unordered_map>
#include <string>
#include <algorithm>


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
@return ret extracted string.
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
int findIndex(std::string query, std::vector<std::string> v) {
	auto index = std::find(v.begin(), v.end(), query);
	if (index != v.end())
		return -1;
	else
		return index - v.begin();
}

/**
Extracts string from a MemoryMapped class

@param charArray MemoryMapped to extract string from.
@param start Index to start extracting string from.
@param end Index to stop extracting string from.
@return ret extracted string.
*/
std::vector<std::string> parseHeader(MemoryMapped &vcf, int &pos) {
	int lastPos = 0;

	//TODO: possible point of failure if file is not formatted properly
	//look for last header line
	while (true) {

		if (vcf[pos] == '\n') {
			if (vcf[pos + 1] != '#')
				break;
			lastPos = pos + 2;
		}
		pos++;
	}

	pos = lastPos;
	std::vector<std::string> colNames;

	//parse the header
	while (true) {
		if (vcf[pos] == '\t') {
			colNames.push_back(getString(vcf, lastPos, pos));
			lastPos = pos + 1;
		}
		else if (vcf[pos] == '\n') {
			colNames.push_back(getString(vcf, lastPos, pos));
			pos++;
			break;
		}
		pos++;
	}

	return colNames;
}


/**
Seperates the case and control IDs

@param vcfDir Full directory path of the VCF file.
@param vcfDir Full directory path of the caseID file.
@param ncolID The number of columns before the sample IDs start in the last line of the headers in the VCF file.
@return A map with sample names as keys and true (case) or false (control) as values.
*/
std::unordered_map<std::string, bool> getIDs(std::string vcfDir, std::string caseIDDir, int ncolID) {

	//open VCF file and find the line with column names
	MemoryMapped vcf(vcfDir);
	int pos = 0;
	std::vector<std::string> IDs = parseHeader(vcf, pos);

	int lineCounter = 0;

	//open caseID file
	MemoryMapped caseID(caseIDDir);
	pos = 0;
	int startPos = pos;
	std::vector<std::string> caseIDs;

	//extract every ID from caseID file, assumes one ID per line
	std::string ID;

	while (pos < caseID.mappedSize()) {
		if (caseID[pos] == '\n') {

			ID = trim(getString(caseID, startPos, pos));
			if (ID.size() > 0) {
				caseIDs.push_back(ID);
				lineCounter++;
				startPos = pos + 1;
			}
		}
		pos++;
	}

	//last line may not end with return character
	std::string lastLine = trim(getString(caseID, startPos, pos));
	if (lastLine.size() > 0) {
		caseIDs.push_back(lastLine);
		lineCounter++;
	}

	//map each ID to true or false:
	//control -> false
	//case -> true
	std::unordered_map<std::string, bool> IDmap;
	int countControl = 0;

	for (size_t i = 0; i < IDs.size(); i++) {

		int index = findIndex(IDs[i], caseIDs);
		if (index > -1) {
			IDmap.insert(std::pair<std::string, bool>(IDs[i], false));
			countControl++;
		}
		else {
			IDmap.insert(std::pair<std::string, bool>(IDs[i], true));
		}
	}

	int countCase = IDmap.size() - countControl;

	std::cout << "This VCF includes " + std::to_string(IDmap.size()) + " samples with " + std::to_string(countControl) +
		" controls and " + std::to_string(countCase) + " cases.\n";

	if (lineCounter != countCase)
		std::cout << "Warning: " + std::to_string(lineCounter) + " lines were counted in case ID file but only " +
		std::to_string(countCase) + " were found to correspond to columns in .vcf file.\n";

	return IDmap;
}