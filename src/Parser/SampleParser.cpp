#include "Parser.h"

static const std::string ERROR_SOURCE = "SAMPLE_INFO_PARSER";

/**
Parses a data file containing tab-separated info for each sample. Verifies the ID column
 is consistent with data from the VCF file.

@param sampleDir Directory of tab-separated data file.
@param IDmap Map that gives a unique int for each sample name.

@return bool if IDs are valid
*/
bool validateSampleIDs(std::string sampleDir, std::map<std::string, int> &IDmap){

    File sampleInfo;
    sampleInfo.open(sampleDir);

    std::vector<std::string> lineSplit;

    std::map<std::string, int> ID;

    while (sampleInfo.hasNext()) {

        lineSplit = split(sampleInfo.nextLine(), SAMPLE_SEP);

        if (lineSplit.size() < 2)
            break;

        std::string sampleID = lineSplit[ID_COL];

        if(ID.count(sampleID) > 0){
            std::string message = "Line " + std::to_string(sampleInfo.lineNumber) +
            " in sample information file - ID not unique in file.";
            throwError(ERROR_SOURCE, message, sampleID);
        }
        if (IDmap.count(sampleID) < 1) {
            std::string message = "Line " + std::to_string(sampleInfo.lineNumber) +
            " in sample information file - ID does not correspond to an ID in VCF file.";
            throwError(ERROR_SOURCE, message, sampleID);
        }

        ID[sampleID] = 1;
    }

    sampleInfo.close();

    std::vector<std::string> missingSample;
    for (auto it = IDmap.begin(); it != IDmap.end(); ++it) {
        std::string id = it->first;

        if (ID[id] != 1)
            missingSample.push_back(id);
    }

    if (missingSample.size() > 0) {

        std::string message = "Samples were found in the provided VCF but are missing in the sample information file:";
        for (int i = 0; i < missingSample.size(); i++)
            message += " '" + missingSample[i] + "' ";

        throwError(ERROR_SOURCE, message);
    }

    return true;
}

/**
Parses a data file containing tab-separated info for each sample. Extracts phenotype column as vector.

@param sampleDir Directory of tab-separated data file.
@param IDmap Map that gives a unique int for each sample name.

@return Vector of phenotype data
*/
VectorXd parseSamplePhenotype(std::string sampleDir, std::map<std::string, int> &IDmap){

    File sampleInfo;
    sampleInfo.open(sampleDir);

    VectorXd Y = VectorXd(IDmap.size());
    for(int i = 0; i<Y.rows(); i++)
        Y[i] = NAN;

    std::vector<std::string> lineSplit;

    while (sampleInfo.hasNext()) {

        lineSplit = split(sampleInfo.nextLine(), SAMPLE_SEP);

        if (lineSplit.size() < 2)
            break;

        std::string sampleID = lineSplit[ID_COL];
        int index = IDmap[sampleID];

        try {
            Y[index] = std::stod(lineSplit[PHENOTYPE_COL]);
        }
        catch (...) {
            if(!lineSplit[PHENOTYPE_COL]=="NA"){
                std::string message = "Line " + std::to_string(sampleInfo.lineNumber) +
                " in sample information file - Unexpected value non-numeric value in phenotype column. Use NA if missing.";
                sampleInfo.close();
                throwError(ERROR_SOURCE, message, lineSplit[PHENOTYPE_COL]);
            }
        }
    }

    sampleInfo.close();
    return Y;
}

/**
Parses a data file containing tab-separated info for each sample. Extracts group IDs as vector.

@param sampleDir Directory of tab-separated data file.
@param IDmap Map that gives a unique int for each sample name.

@return Vector of group ID
*/
VectorXd parseSampleGroupID(std::string sampleDir, std::map<std::string, int> &IDmap){

    File sampleInfo;
    sampleInfo.open(sampleDir);

    std::map<std::string, int> groupIDMap;
    int groupIndex = 0;
    VectorXd G = VectorXd(IDmap.size());
    for(int i = 0; i < G.rows(); i++)
        G[i] = -1;

    std::vector<std::string> lineSplit;
    while (sampleInfo.hasNext()) {

        lineSplit = split(sampleInfo.nextLine(), SAMPLE_SEP);

        if (lineSplit.size() < 2)
            break;

        std::string sampleID = lineSplit[ID_COL];
        std::string groupID = lineSplit[GROUP_COL];
        int index = IDmap[sampleID];

        if (!groupIDMap.count(groupID)){
            groupIDMap[groupID] = groupIndex;
            groupIndex++;
        }

        G[index] = groupIDMap[groupID];
    }

    sampleInfo.close();
    return G;
}

/**
Parses a data file containing tab-separated info for each sample. Extracts read depth status as map.

@param sampleDir Directory of tab-separated data file.
@param IDmap Map that gives a unique int for each sample name.
@param highLowCutOff Threshold between high and low read depth.
@param G Vector of group IDs.

@return Map from group ID to read depth
*/
std::map<int, Depth> parseSampleReadDepth(std::string sampleDir, std::map<std::string, int> &IDmap, VectorXd &G, int highLowCutOff){

    File sampleInfo;
    sampleInfo.open(sampleDir);

    std::map<int, Depth> groupDepth;

    std::vector<std::string> lineSplit;

    while (sampleInfo.hasNext()) {

        lineSplit = split(sampleInfo.nextLine(), SAMPLE_SEP);

        if (lineSplit.size() < 2)
            break;

        std::string sampleID = lineSplit[ID_COL];
        int index = IDmap[sampleID];
        int groupIndex = G[index];

        if (!groupDepth.count(groupIndex)) {

            std::string depth = lineSplit[DEPTH_COL];

            if (depth.find("H") != std::string::npos ||
                depth.find("h") != std::string::npos)
                groupDepth[groupIndex] = Depth::HIGH;
            else if (depth.find("L") != std::string::npos ||
                depth.find("l") != std::string::npos)
                groupDepth[groupIndex] = Depth::LOW;
            else {
                try {
                    int d = std::stoi(depth);
                    if(d >= highLowCutOff)
                        groupDepth[groupIndex] = Depth::HIGH;
                    else
                        groupDepth[groupIndex] = Depth::LOW;
                }
                catch (...) {
                    std::string message = "Line " + std::to_string(sampleInfo.lineNumber);
                    message += " in sample information file: Unexpected value in read depth column. Should be a numeric value or one of 'L', 'H'.";
                    sampleInfo.close();
                    throwError(ERROR_SOURCE, message, lineSplit[DEPTH_COL]);
                }
            }
        }
    }

    sampleInfo.close();
    return groupDepth;
}


/**
Parses a data file containing tab-separated info for each sample.

@param covariate Column of covariates read from sample data file.
@param isNumeric True if covariate is a numeric value and not categorical.
@param columnIndex Index of the covariate column.
@return Matrix of covariate values
*/
MatrixXd handleCovariate(std::vector<std::string> covariate, bool isNumeric, int columnIndex) {

    if (isNumeric) {
        VectorXd z(covariate.size());
        for (int i = 0; i < covariate.size(); i++) {

            if (covariate[i] == "NA" || covariate[i] == "")
                z[i] = NAN;
            else
                z[i] = stod(covariate[i]);
        }
        return z;
    }
    else {
        std::map<std::string, int> covIDMap;
        int id = 0;

        for (int i = 0; i < covariate.size(); i++) {
            if (!covIDMap.count(covariate[i])) {
                covIDMap[covariate[i]] = id;
                id++;
            }
        }

        if (covIDMap.size() > 5) {
            std::string message = "In sample data file: Covariate in column ";
            message += std::to_string(columnIndex + 5) + " is being treated as categorical and has ";
            message += std::to_string(covIDMap.size()) + " unique values. This could have a negative impact on computation time.";
            printWarning(message);
        }

        MatrixXd z(covariate.size(), covIDMap.size() - 1);
        z = z.setZero();

        int col;
        for (int i = 0; i < covariate.size(); i++) {

            if (covariate[i] == "NA" || covariate[i] == "") {
                for (int j = 0; j < covIDMap.size() - 1; j++)
                    z(i, j) = NAN;
            }
            else {
                col = covIDMap[covariate[i]];

                //don't include last categorical covariate
                if(col < id-1)
                    z(i, col) = 1;
            }
        }
        return z;
    }
}

/**
Parses a data file containing tab-separated info for each sample.

@param covariates Column(s) of covariates read from sample data file.
@return Matrix of covariate values
*/
MatrixXd handleCovariates(std::vector<std::vector<std::string>> covariates) {
    MatrixXd Z;

    for (int i = 0; i < covariates[0].size(); i++) {

        std::vector<std::string> covariate;
        bool isNumeric = true;

        for (int j = 0; j < covariates.size(); j++) {
            covariate.push_back(covariates[j][i]);
            if (isNumeric) {
                try {
                    std::stod(covariates[j][i]);
                }
                catch (...) {
                    std::string cov = covariates[j][i];
                    if (cov != "NA" && cov != "")
                        isNumeric = false;
                }
            }
        }
        MatrixXd cov = handleCovariate(covariate, isNumeric, i);

        if (i == 0)
            Z = cov;
        else {
            MatrixXd temp = Z;
            MatrixXd T(temp.rows(), temp.cols() + cov.cols());
            Z=T;
            Z << temp, cov;
        }
    }

    MatrixXd interceptZ(Z.rows(), Z.cols()+1);
    interceptZ << Eigen::VectorXd::Constant(Z.rows(), 1), Z;
    return interceptZ;
}

/**
Parses a data file containing tab-separated info for each sample. Extracts covariates as a matrix.

@param sampleDir Directory of tab-separated data file.
@param IDmap Map that gives a unique int for each sample name.

@return Matrix of covariates of group ID
*/
MatrixXd parseSampleCovariates(std::string sampleDir, std::map<std::string, int> &IDmap){

	File sampleInfo;
    sampleInfo.open(sampleDir);
	    
    std::vector<std::vector<std::string>> covariates(IDmap.size());

	std::vector<std::string> lineSplit;

	while (sampleInfo.hasNext()) {

        lineSplit = split(sampleInfo.nextLine(), SAMPLE_SEP);

		if (lineSplit.size() < 2)
			break;

        std::string sampleID = lineSplit[ID_COL];
        int index = IDmap[sampleID];

		std::vector<std::string> cov;
        for (size_t i = COV_COL; i < lineSplit.size(); i++)
			cov.push_back(lineSplit[i]);

        covariates[index] = cov;
	}

    MatrixXd Z;
	if(covariates.size() > 0)
        Z = handleCovariates(covariates);

	sampleInfo.close();

    return Z;
}