#pragma once
#include "stdafx.h"
#include "MemoryMapped.h"

#include <string>
#include <vector>
#include <unordered_map>


std::unordered_map<std::string, bool> getIDs(std::string vcfDir, std::string caseIDDir, int ncolID);
inline std::string trim(std::string str);
std::vector<std::string> parseHeader(MemoryMapped &vcf, int &pos);
inline std::string getString(MemoryMapped &charArray, int start, int end);
int findIndex(std::string query, std::vector<std::string> v);