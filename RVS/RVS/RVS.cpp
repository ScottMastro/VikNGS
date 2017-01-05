#include "stdafx.h"
#include "MemoryMapped.h"

#include <iostream>  
#include <string>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <iomanip>


using namespace std;

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
	string chr;
	int loc;

	vector<genotypeLikelihood> gl;
	vector<double> p;
	vector<double> eg;
	double maf = NULL;

	inline bool operator<(SNP& snp2) { return this->loc < snp2.loc; }
	inline double L00(int index) { return this->gl[index].L00; }
	inline double L01(int index) { return this->gl[index].L01; }
	inline double L11(int index) { return this->gl[index].L11; }

};

inline bool locCompare(SNP lhs, SNP rhs) { return lhs < rhs; }

//========================================================


//Extracts string from a MemoryMapped class
inline string getString(MemoryMapped &charArray, int start, int end) {
	string ret;
	for (; start < end; start++) { ret += charArray[start]; }
	return ret;
}

string trim(string str)
{
	str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
	return str;
}

int findIndex(string query, vector<string> v) {
	auto index = std::find(v.begin(), v.end(), query);
	if (index != v.end())
		return -1;
	else 
		return index - v.begin();
}

vector<string> parseHeader(MemoryMapped &vcf, int &pos) {
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
	vector<string> colNames;

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


/* 
#' Function to seperates the case and control IDs.
#'
#' This function read the header of the vcf file (rows starting with # in vcf) and split the case and controls according to the caseID_file. It returns two vectors each of which is the position index of case or control in the list of samples in vcf file (the last row starting with #).
#' @param vcf_file The saving address and filename  of the vcf file in the format of 'address/filename'.
#' @param caseID_file The  file including case IDs in the format of 'address/filename'.
#' @param nhead The number of the lines in the vcf file which started with #.
#' @param ncol_ID The number of columns before the sample IDs start in the last line of the headers in the vcf file, default=9.
#' @return a list which inclues two vectors that correspond to the position index of the cases or controls in the last head row of the VCF with IDs, but columns are shifted by the ncol_ID columns.
# @examples getIDs("../data/1g113low_1g56exomehigh_filtered.hg19.chr11_1000snps.vcf","../data/bam.txt",128,9)
*/
std::unordered_map<string, bool> getIDs(string vcfDir, string caseIDDir, int ncolID) {

	MemoryMapped vcf(vcfDir);
	int pos = 0;
	std::vector<string> IDs = parseHeader(vcf, pos);

	int lineCounter = 0;

	MemoryMapped caseID(caseIDDir);
	pos = 0;
	int startPos = pos;
	std::vector<string> caseIDs;

	while (pos < caseID.mappedSize()) {
		if (caseID[pos] == '\n') {
			caseIDs.push_back(trim(getString(caseID, startPos, pos)));
			lineCounter++;
			startPos = pos + 1;
		}
		pos++;
	}
	string lastLine = trim(getString(caseID, startPos, pos));
	if (lastLine.size() > 0) {
		caseIDs.push_back(lastLine);
		lineCounter++;
	}

	std::unordered_map<string, bool> IDmap;

	int countControl = 0;

	//control -> false
	//case -> true
	for (size_t i = 0; i < IDs.size(); i++) {

		int index = findIndex(IDs[i], caseIDs);
		if (index > -1) {
			IDmap.insert(std::pair<string, bool>(IDs[i], false));
			countControl++;
		}
		else {
			IDmap.insert(std::pair<string, bool>(IDs[i], true));
		}
	}

	int countCase = IDmap.size() - countControl;

	cout << "This VCF includes "  + to_string(IDmap.size()) + " samples with " + to_string(countControl) +
		" controls and " + to_string(countCase) + " cases.\n";

	if (lineCounter != countCase)
		cout << "Warning: " + to_string(lineCounter) + " lines were counted in case ID file but only " +
		to_string(countCase) + " were found to correspond to columns in .vcf file.\n";

	return IDmap;
}

/*
#' Function to get Likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields
#'
#' This function fetches the genotype likelihoods from the VCF file (wo header) for the analysis.
#' The first step is to read the vcf file as a matrix and outputs an R object with 6 lists including 5 vectors and 1 list including 3 matrix:
#' 5 vectors: chr (1st column of VCF file), position (2nd column of VCF file, POS),ref allele (4th column, REF), alt allele (5th column, ALT), filter indicator (7th column, FILT)
#' 1 list: likelihood of L(00|ref,alt), L(01|ref,alt), L(11|ref,alt), which are three matrix with row for sample and col for variants.
#' They are calculated by \code{\link{get_lsingle_vcf}} from the PL part of the 10th column and onwards of VCF in format of GT:AD:DP:GQ:PL, example: 0/0:2,0:4:6:0,6,42
#' GL = genotype likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields in the form of log10-scale,
#' also called as Phred-scaled likelihoods: -log10(L(00|ref,alt)), -log10(L(01|ref,alt)), -log10(L(11|ref,alt)).
#' PL : the phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field) (Integers)
#'
#' @param M  The vcf file without the headers in a matrix, each row stands for a variant with number of columns = ncol_ID + number of sampls.
#' @param ncol_ID the number of columns before the data for each individual
#' @return Outputs an R object with
#' 1. vector for chromosome of variants
#' 2. vector for the position of variant location
#' 3. vector for the reference allele
#' 4. ALT allele (vector)\cr
#' 5. FILTER indicator for 'PASS' (vector)\cr
#' 6. likelihoods: L(D|0), L(D|1), L(D|2). List of 3 matrix:
#' L0 - matrix of P(D|0), rows are individuals and columns are SNPS\cr
#' L1 - matrix of P(D|1), rows are individuals and columns are SNPS\cr
#' L2 - matrix of P(D|2), rows are individuals and columns are SNPS\cr
*/
SNP buildSNP(MemoryMapped &vcf, vector<int> ind, int ncolID) {

	SNP snp;
	snp.chr = getString(vcf, ind[0], ind[1] - 1);
	snp.loc = stoi(getString(vcf, ind[1], ind[2] - 1));

	int indexPL = 0;
	int pos = ind[8];
	while (true) {
		if (vcf[pos] == ':')
			indexPL++;
		else if (vcf[pos] == 'P' && vcf[pos + 1] == 'L')
			break;
		else if (vcf[pos] == '\t') {
			cout << "ERROR: Phred-scaled likelihoods (PL) not found in FORMAT column ";
			cout << "or FORMAT is not the 9th column in the .vcf file.";
		}
		pos++;
	}

	vector<genotypeLikelihood> likelihoods;
	int startPos;
	int counter;
	double l00;
	double l01;
	double l11;

	for (int i = ncolID; i < ind.size(); i++) {
		pos = ind[i];
		counter = 0;

		if (vcf[pos] == '.') {
			genotypeLikelihood gl;
			gl.L00 = NULL;
			gl.L01 = NULL;
			gl.L11 = NULL;
			likelihoods.push_back(gl);

			continue;
		}

		while (true) {
			if (vcf[pos] == ':') {
				counter++;
				if (counter == indexPL) {
					pos++;
					break;
				}
			}
			pos++;
		}
		
		startPos = pos;

		while (vcf[pos] != ',')
			pos++;
		pos++;
		l00 = stod(getString(vcf, startPos, pos));
		startPos = pos;

		while (vcf[pos] != ',')
			pos++;
		pos++;
		l01 = stod(getString(vcf, startPos, pos));
		startPos = pos;

		while (vcf[pos] != '\t' && vcf[pos] != ':' && vcf[pos] != '\n')
			pos++;
		l11 = stod(getString(vcf, startPos, pos));

		genotypeLikelihood gl;
		gl.L00 = pow(10, -l00*0.1);
		gl.L01 = pow(10, -l01*0.1);
		gl.L11 = pow(10, -l11*0.1);
		likelihoods.push_back(gl);
	}

	snp.gl = likelihoods;

	return snp;
}

vector<SNP> filterVCFSNPs(string vcfDir, int ncolID, double missingTh, std::unordered_map<string, bool> &IDmap) {

	MemoryMapped vcf(vcfDir);
	vector<SNP> snps;

	int pos = 0;
	vector<string> colNames = parseHeader(vcf, pos);
	size_t ncol = colNames.size();

	//TODO: point of failure: columns not ordered properly
	int chrom = 0;
	int loc = 1;
	int ref = 3;
	int alt = 4;
	int filter = 6;


	int ncase = 0;
	int ncontrol = 0;
	for (auto& p : IDmap) {
		ncase += p.second;
		ncontrol += !p.second;
	}

	int allowedMissCase = floor(ncase * missingTh);
	int allowedMissControl = floor(ncontrol * missingTh);

	int filterPass = 0;
	int filterBase = 0;
	int filterCase = 0;
	int filterControl = 0;
	int filterDuplicate = 0;

	//can parallelize if desired:
	//======================================
	int start = pos;
	int end = vcf.mappedSize();
	
	int currInd = 0;
	int missCase = 0;
	int missControl = 0;
	bool filtered = false;

	int counter = 0;
	int successcounter = 0;
	while (start < end) {
		vector<int> ind;
		ind.push_back(start);
		currInd = 1;
		missCase = 0;
		missControl = 0;
		filtered = false;

		while (true) {
			//parsing the columns before ID columns
			if (currInd <= ncolID) {
				if (vcf[start] == '\t') {
					ind.push_back(start + 1);

					//filter out SNPs where ref and alt are not single base
					if (currInd == 3 || currInd == 4) {
						if (vcf[start + 2] != '\t') {
							filtered = true;
							filterBase++;
							break;
						}
					}

					//filter out SNPs where filter is not "PASS"
					else if (currInd == 6) {
						if (vcf[start + 1] != 'P' ||
							vcf[start + 2] != 'A' ||
							vcf[start + 3] != 'S' ||
							vcf[start + 4] != 'S') {
							filtered = true;
							filterPass++;
							break;
						}
					}
					currInd++;
				}
			}
			//parsing ID columns
			else {
				if (start < vcf.mappedSize() && vcf[start] == '\t') {
					ind.push_back(start + 1);
					if (vcf[start + 1] == '.') {
						if (IDmap.at(colNames[currInd])) {
							missCase++;
							if (missCase > allowedMissCase) {
								filtered = true; 
								filterCase++;
								break;
							}
						}
						else {
							missControl++;
							if (missControl > allowedMissControl) {
								filtered = true;
								filterCase++;
								break;
							}
						}
					}
					currInd++;
				}

				else if (start >= vcf.mappedSize() || vcf[start] == '\n')
					break;
			}
			start++;
		}

		while (start - 1 < vcf.mappedSize() && vcf[start - 1] != '\n')
			start++;

		counter++;
		if (!filtered) {
			successcounter++;
			snps.push_back(buildSNP(vcf, ind, ncolID));

		}
	}
	//======================================

	std::sort(snps.begin(), snps.end(), locCompare);

	set<size_t> toDelete;
	size_t j;

	for (size_t i = 0; i < snps.size(); i++) {
		j = 1;
		while (i + j < snps.size() && snps[i].loc == snps[i+j].loc) {

			if (snps[i].chr == snps[i + j].chr) {
				toDelete.insert(i);
				toDelete.insert(i + j);
			}
			j++;
		}
	}

	for (auto it = toDelete.rbegin(); it != toDelete.rend(); it++) {
		snps.erase(snps.begin() + *it);
		filterDuplicate++;
	}

	cout << filterPass;
	cout << " SNPs were filtered by FILTER column\n";
	cout << filterBase;
	cout << " SNPs were filtered by ALT and REF column\n";
	cout << filterCase;
	cout << " SNPs were filtered by missing case values\n";
	cout << filterControl;
	cout << " SNPs were filtered by missing control values\n";
	cout << filterDuplicate;
	cout << " SNPs were filtered by duplication\n";
	cout << successcounter;
	cout << " SNPs remain after filtering\n";

	return snps;
}

/*
#' Function to use EM algorithm to estimate the genotype frequencies in the sample
#'
#' Given a matrix (M) with dimension number of sample by 3 consisting of genotype likelihoods for each individual, this function uses the EM algorithm to estimate genotype frequencies.
#'  It returns a vector with three decimals to stands for probability of 0, 1 or 2 minor alleles in the whole sample
#'
#' @param M  matrix of genotype likelihoods derived from Phred-scale likelihood in VCF input for one variant.
#' @return a vector with three decimals to stands for probability of 0, 1 or 2 minor alleles.
#' @export
*/
vector<double> calcEM(SNP &snp) {
	double p = 0.15;
	double q = 0.15;
	double qn = 1;
	double pn = 0;
	double dn = 0;
	double d = 0;

	int glCounter=0;
	
	double Ep;
	double Eq;
	double pD;

	int k = 0;
	

	while (pow((pn - p), 2) + pow((qn - q),2) > 0.000001) {

		d = 1 - p - q;

		glCounter = 0;
		Ep = 0;
		Eq = 0;


		for (size_t i = 0; i < snp.gl.size(); i++) {
			//i = snp.gl.size() - 1 - i;
			if (snp.L00(i) == NULL)
				continue;

			glCounter++;

			pD = 1/(p * snp.L00(i) + q * snp.L01(i) + d * snp.L11(i));
			Ep += p * snp.L00(i) * pD;
			Eq += q * snp.L01(i) * pD;
		}

		pn = p;
		qn = q;
		dn = 1 - q - p;
		p = Ep / glCounter;
		q = Eq / glCounter;
		
		k++;
		if (k == 1000)
			break;
	}

	vector<double> freq;
	freq.push_back(p);
	freq.push_back(q);
	freq.push_back(1 - p - q);

	return freq;
}


/*
#' Function to calculate the conditional expected  genotype probability \eqn{E(P(G_{ij}|D_{ij}))} given the genotype likelihoods \eqn{P(D_{ij}|G_{ij}=g)} and frequencies. Simplified version of \code{\link{calc_EG_general}}.
#'
#' Calculate the conditional expected genotype probability given the observed sequencing data \eqn{D_{ij}} from formula \eqn{E(P(G_{ij}|D_{ij}))=\sum_{g=0}^2 gP(G_{ij}=g|D_{ij})}, where \eqn{ P(G_{ij}=g|D_{ij})=P(D_{ij}|G_{ij}=g)*P(G_{ij}=g)/P(D_{ij})}.
#' \eqn{P(D_{ij}|G_{ij}=g)} (input M) is from the VCF file or function \code{\link{calc_pobs_ndepth}}, \eqn{p(G_{ij}=g)} (input p) is output from function \code{\link{calc_EM}}.
#' All values for \eqn{P(G_{ij}=g|D_{ij})} are scaled by \eqn{P(D_{ij})}. 
## It is called by \code{\link{generate_seqdata}}.
## talk to Andriy, he said in the kk loop, he did it for more general case, we only need rdv=1 for our case.
## so I renamed Andriy's get_EG to calc_EG_general, delete the kk loop in calc_EG function.
## all values are rescaled by P(D_ij), (see m below, but it is not probability, so m/sum(m) gives us sth. like prob.
#'
#' @param M   genotype likelihoods with dimension of number of sample times 3 (double), each column is the genotype likelihood \eqn{P(D_{ij}|G=AA, Aa\ or\ aa)} for one locus. It uses output from \code{\link{get_exp_geno}} or \code{\link{get_exp_MAF}} for VCF input.
#' @param p   genotype frequencies \eqn{P(G=AA, Aa\ or\ aa)} or \eqn{p(G=0,1,2\ minor\ allele)} for each SNP, it is a vector of length 3 and uses the output from  function \code{\link{calc_EM}}.
#' @param rdv   read depth for all samples. Dummy variable in \code{\link{calc_EG}}, listed here to be consitent with \code{\link{calc_EG_general}}. It is a vector of integers of one with length equal to the number of samples.
#' @return  a vector with the same length as rdv, containing conditional expectation probability \eqn{ E(P(G_{ij}|D_{ij}))}.
*/
vector<double> calcEG (SNP &snp, vector<double> &p) {
	
	double m0;
	double m1;
	double m2;
	double m;

	vector<double> EG;

	for (size_t i = 0; i < snp.gl.size(); i++) {

		if (snp.L00(i) != NULL) {
			m0 = snp.L00(i) * p[0];
			m1 = snp.L01(i) * p[1];
			m2 = snp.L11(i) * p[2];
			m = 1 / (m0 + m1 + m2);

			EG.push_back(m1*m + 2 * m2*m);
		}
		else
			EG.push_back(NULL);

	}

	return EG;
}

/*
#' Function to generates the expected probabilities of the genotypes E(G_ij|D_ij).
#'
#' Using the genotype likelihood for case and controls from VCF file to generates the population frequency by calling function \code{\code{calc_EM}} and then use it to calculate the expected genotype probabilities E(G_ij|D_ij) by calling function \code{\link{calc_EG}}. Variants with homozygous call in the whole sample (the std of their E(G_ij|D_ij) <10^4) will be removed and the expected genotype probability will be reorganized so that the top rows are for case and the rest for controls.
#' @param A0M Genotype likelihood for major homozygous of the cases, output of filter_VCF_SNPs.
#' @param A1M Genotype likelihood for minor heterozygous of the cases, output of filter_VCF_SNPs.
#' @param A2M Genotype likelihood for minor homozygous of the cases, output of filter_VCF_SNPs.
#' @param B0M Genotype likelihood for major homozygous of the controls, output of filter_VCF_SNPs.
#' @param B1M Genotype likelihood for minor heterozygous of the controls, output of filter_VCF_SNPs.
#' @param B2M Genotype likelihood for minor homozygous of the controls, output of filter_VCF_SNPs.
#' @param chr chromosome of the variant, output from the  function filter_VCF_SNPs
#' @param loc location of the variant, output from the function filter_VCF_SNPs.
#' @return a list includes a matrix (exp_cond_prob) contains the expected genotype likelihood for each person (row) at each locus (column); a matrix (pop_freq) for genotype frequency in all the sample at each locus (size: number of variant *3); a data frame (SNPs) contains chr and loc of the returned variants and number of variants (nrm_hom) removed because of homozygous calls in the whole sample.
*/
void getExpGeno(vector<SNP> &snps) {

	vector<double> EG;
	vector<size_t> filterIndex;

	double mean;
	double sdSum;

	double n;
	int filterCount = 0;

	for (size_t i = 0; i < snps.size(); i++) {
		vector<double> p;

		p = calcEM(snps[i]);
		snps[i].p = p;
		EG = calcEG(snps[i], p);
		snps[i].eg = EG;

		//check and filter if variant is homozygous
		mean = 0;
		n = 0;
		for (size_t j = 0; j < EG.size(); j++) {
			if (EG[j] != NULL) {
				mean += EG[j];
				n++;
			}
		}
		mean = mean / n;

		sdSum = 0;
		for (size_t j = 0; j < EG.size(); j++)
			if (EG[j] != NULL)
				sdSum += pow((EG[j] - mean), 2);

		if (1e-8*(n - 1) > sdSum)
			filterIndex.push_back(i);
	}

	for (size_t i = 0; i < filterIndex.size(); i++) {
		snps.erase(snps.begin() + filterIndex[i]);
		filterCount++;
	}

	if (filterCount > 0) {
		cout << filterCount;
		cout << " SNPs were removed because of homozygous call in all samples\n";
	}
	return;
}


void getExpMAF(vector<SNP> &snps, double mafCut, bool common) {
	double maf;
	int mafCounter = 0;


	for (size_t i = 0; i < snps.size(); i++) {
		maf = 0.5 * snps[i].p[1] + snps[i].p[2];
		if (maf > 0.5)
			maf = 1 - maf;

		if (common) {
			if (maf >= mafCut) {
				mafCounter++;
				snps[i].maf = maf;
			}
		}
		else {
			if (maf < mafCut) {
				mafCounter++;
				snps[i].maf = maf;
			}
		}
	}

	cout << mafCounter;
	cout << " out of ";
	cout << snps.size();
	cout << " variants satisfy the MAF condition provided\n";
}

/*
#' the main function to call the VCF file and return the  expected probabilities of the genotypes E(G_ij|D_ij) from VCF file
#'
#' @param file the VCF file
#' @param caseID_file the caseID files
#' @param nhead the row number of your VCF which contains all sample ID
#' @param ncol_ID the number of columns before the sample ID
#' @param missing_th the missing rate cut off for a SNP to be filtered, if missing>missing_th in case or controls, the SNP will be removed
#' @param nvar_tot should be the number of SNPs in your VCF file, but if VCF is too large, you can choose only read in the first nsnp SNPs
#' @param nread the number of SNPs read each time by R (nread<=nsnp)
#' @param maf_cut the MAF to be used to filter out SNP,
#' @param common Boolean variable to indicate common (TRUE) or rare (FALSE) variants 
#'
#' @return a list  for case or control (depends on group) including expected conditional probabilities, SNP info, Y, ncase/ncont etc 
#' @export
*/
void vcf_process(string file, string caseID_file) {

	std::unordered_map<string, bool> IDmap = getIDs(file, caseID_file, 9);
	std::vector<SNP> snps = filterVCFSNPs(file, 9, 0.2, IDmap);

	if (snps.size() == 0)
		cout << "No SNPs left after filtering .vcf file\n";


	/*//print out genotype likelihood

	for (size_t i = 0; i < snps[45].gl.size(); i++) {
		cout << snps[45].L00(i);
		cout << '\t';
		cout << snps[45].L01(i);
		cout << '\t';
		cout << snps[45].L11(i);
		cout << '\n';
	}*/

	getExpGeno(snps);
	
	if (snps.size() == 0)
		cout << "No SNPs left removing homozygous variants\n";


	bool common = true;
	double mafCut = 0.05;

	getExpMAF(snps, mafCut, true);


	cout << "Chr\tLoc\tMAF\n";

	for (size_t i = 0; i < snps.size(); i++) {

		if (snps[i].maf != NULL) {
			cout << snps[i].chr;
			cout << '\t';
			cout << snps[i].loc;
			cout << '\t';
			cout << snps[i].maf;
			cout << '\n';

		}
	}


	return;
}


/*
vcf_process<-function(file = 'example/example_1000snps.vcf',
	caseID_file = 'example/caseID.txt',
	nhead = 128, ncol_ID = 9, missing_th = 0.2, nvar_tot = 1000, nread = 300, maf_cut = 0.05, common = TRUE)
{
	geno<-get_exp_MAF(exp_prob$pop_frq, exp_prob$exp_cond_prob, exp_prob$SNPs[, 'chr'], exp_prob$SNPs[, 'loc'], ncase = length(IDs$cases), maf_cut = maf_cut, common = common)

		if (!class(geno) % in%'NULL') {
			ngeno = nrow(geno$SNPs)
				if (common == TRUE) {
					group = 'common'
				}
				else { group = 'rare' }
				cat('There are ', ngeno, group, 'variants. \n\n')

					SNPs = geno$SNPs; P = geno$P; Geno = geno$Geno; ncase = geno$ncase
					ncont = geno$ncont; nsnp = geno$nsnp; Y = geno$Y;
				write.table('SNPs', paste0(nsnp, 'snps_left_for', ncase, 'ncase', ncont, 'ncont.txt', row.names = F, quote = F, sep = '\t', colnames = T))
					return(geno)
		}


}
*/


int main() {

	vcf_process("C:/Users/Scott/Desktop/RVS-master/example/example_1000snps.vcf",
				"C:/Users/Scott/Desktop/RVS-master/example/caseID.txt");



	while (true) {

	}

	return 0;

}
