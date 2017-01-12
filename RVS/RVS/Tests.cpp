#include "stdafx.h"
#include "RVS.h"

#include <iostream>
#include <vector>


/*
Calculates the robust variance of E(G | D). var(x) = E(x^2) - E(x)^2

@param p Genotype frequency from EM algorithm
@return Robust variance of E(G | D)
*/
inline double calcRobustVar(std::vector<double> p) {
	return (4 * p[2] + p[1]) - pow(2 * p[2] + p[1], 2);
}


void calcMeanVar(std::vector<bool> &IDmap, std::vector<SNP> &snps) {
	double var;
	double mean;
	double n;
	double controlvar;
	double controlmean;
	double ncontrol;
	double casevar;
	double casemean;
	double ncase;
	double eg;

	for (size_t j = 0; j < snps.size(); j++) {

		var = 0;
		mean = 0;
		n = 0;
		controlvar = 0;
		controlmean = 0;
		ncontrol = 0;
		casevar = 0;
		casemean = 0;
		ncase = 0;

		for (size_t i = 0; i < snps[j].EG.size(); i++) {
			eg = snps[j].EG[i];
			if (eg != NULL) {
				if (!IDmap[i]) {
					ncontrol++;
					controlmean += eg;
				}
				else {
					ncase++;
					casemean += eg;
				}
			}
		}
		
		mean = casemean + controlmean;
		n = ncase + ncontrol;
	
		mean /= n;
		controlmean /= ncontrol;
		casemean /= ncase;

		for (size_t i = 0; i < snps[j].EG.size(); i++) {
			eg = snps[j].EG[i];

			if (eg != NULL) {
				var += pow((eg - mean), 2);

				if (!IDmap[i])
					controlvar += pow((eg - controlmean), 2);
				else
					casevar += pow((eg - casemean), 2);
			}
		}

		var /= n - 1;
		controlvar /= ncontrol - 1;
		casevar /= ncase - 1;

		snps[j].var = var;
		snps[j].mean = mean;
		snps[j].n = n;
		snps[j].controlvar = controlvar;
		snps[j].controlmean = controlmean;
		snps[j].ncontrol = ncontrol;
		snps[j].casevar = controlvar;
		snps[j].casemean = casemean;
		snps[j].ncase = ncase;
	}
	
	return;
}


/*
RVS using asymptotic distribution for score test statistic. This functions includes
two association tests which differ in the estimation of the variance of score. It uses
var_case(E(G | D)) when rvs = false and var_case(G) when rvs = true.

Reference: http://www.ncbi.nlm.nih.gov/pubmed/22570057
Reference: http://www.ncbi.nlm.nih.gov/pubmed/24733292
@param snps Vector of SNPs.
@param IDmap Vector with phenotypes (case/control).
@param rvs Indicates how to calculate variance
@return p-values for the SNPs
*/
std::vector<double> RVSasy(std::vector<SNP> &snps, std::vector<bool> &IDmap, bool rvs) {
	
	std::vector<double> pvals;
	SNP snp;
	double v;
	double p;

	for (size_t j = 0; j < snps.size(); j++) {
		snp = snps[j];

		if (rvs) {
			p = snp.ncase / snp.n;

			v = (1-p) * calcRobustVar(snp.p) + p * snp.controlvar;
			v = snp.var / v;
		}
		else
			v = 1;

		pvals.push_back(chiSquareOneDOF(scoreTest(IDmap, snp.EG) * v));
	}

	return pvals;
}

void RVSbtrap(std::vector<SNP> &snps, std::vector<bool> &IDmap, bool rvs, int nboot) {

	SNP snp = snps[0];

	double maf = 0;
	double casevar;
	double p;

	double s;


	double tObs;


	maf = snp.mean / 2;

	if (maf < 0.05)
		std::cout << "Warning: MAF of the SNP is" + std::to_string(maf) + ", try to group with other rare variants and use RVS_rare1.\n";
	if (1 - maf < 0.05)
		std::cout << "Warning: MAF of the SNP is" + std::to_string(1 - maf) + ", try to group with other rare variants and use RVS_rare1.\n";

	casevar = calcRobustVar(snp.p);

	if (rvs)
		tObs = (snp.casemean - snp.controlmean) / sqrt(casevar / snp.ncase + snp.controlvar / snp.ncontrol);
	else {
		//calculate score test
		p = snp.ncase / snp.n;
		s = (1 - p) * snp.casemean * snp.ncase - p * snp.controlmean * snp.ncontrol;
		tObs = pow(s, 2) / (p * (1 - p) * snp.n * snp.var);
	}
	
	std::cout << tObs;
	
	return;
}


/*

#' Function to calculate the score test for given expected genotype probability for case and controls seperately (M1, M2)
#'
#' This function calcuates the score test \eqn{T=S^2/var(S)}, for a given \eqn{j, S_j=\sum_i (Y_i-bar(Y))E(G_{ij}|D_{ij})=\sum_{case}(1-bar(Y))E(G|D)-\sum_{cont}bar(Y)E(G|D)}. It is called in \code{\link{regScore_perm}}.
#' @param M1   a vector of the expected  genotype probability for case \eqn{E(G_{ij}|D_{ij})} on one snp
#' @param M2   a vector of the expected genotype probability for control \eqn{E(G_{ij}|D_{ij})} on one snp
#' @return the score test statistic calculated from M1 and M2
calc_score_test = function(M1,M2){
X = c(M1,M2)
Y = c(rep(1,length(M1)),rep(0,length(M2)))
p = length(M1)/length(X)
q = 1 - p
S = q*sum(M1)-p*sum(M2)
vs = p*q*(length(X))*var(X)
x = S^2/vs
return (x)
}










#' use RVS to test associaton by bootrtrap, given phenotype (Y), expected values of genotypes for case and controls (X)
#' estimated genotype frequency (P) and number of times of bootstrap (nboot)
#' @param Y - a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X - a vector of genotype for a snp, first ncase and then ncont
#' @param P (vector with length 3  (double)) - estimated genotype frequency for a variant P(G=0), P(G=1) and P(G=2)
#' @param nboot - number of bootstrap
#' @param RVS, logic variable to indicate the estimation of variance of the score statistic. 
#' @return   p-values for the snp
#' @export
RVS_btrap = function(Y, X, P, nboot, RVS = 'TRUE') {
	if (!RVS %in% c('TRUE', 'True', 'true', 'T', 't', 'FALSE', 'False', 'false', 'F', 'f')) {
		cat('Wrong input for option RVS in RVS_brap, should be True or False.\n')
			return(NULL)
	}
	Y = Y[!is.na(X)]
		X = X[!is.na(X)]
		maf = sum(X) / (length(X) * 2); maf0 = min(maf, 1 - maf)
		if (maf0<0.05) { cat('Warning: MAF of the SNP is', maf0, ', try to group with other rare variants and use RVS_rare1.\n') }

	ncase1 = sum(Y == 1)
		ncont1 = sum(Y == 0)
		p = length(X[Y == 1]) / length(X)
		q = 1 - p
#  if(RVS %in% c('TRUE','True','true','T','t')) {
		vcase = as.numeric(calc_robust_Var(P))
#  } else{
		#    vcase = var(X[Y == 1])
#  }
		vcont = var(X[Y == 0])
		if (RVS %in% c('TRUE', 'True', 'true', 'T', 't')) {
			Tobs = (mean(X[Y == 1]) - mean(X[Y == 0])) / sqrt(vcase / ncase1 + vcont / ncont1)
		}
		else {
			X1 = X[Y == 1]; X2 = X[Y == 0]
				Tobs = calc_score_test(X1, X2)
		}
		X1 = X[Y == 1] - mean(X[Y == 1])
			X2 = X[Y == 0] - mean(X[Y == 0])
			C = NULL
			if (RVS %in% c('TRUE', 'True', 'true', 'T', 't')) {
				for (j in 1 : nboot) {
					Xca = sample(X1[], ncase1, replace = TRUE)
						Xco = sample(X2[], ncont1, replace = TRUE)
						vcase = var(Xca)
						vcont = var(Xco)
						C = c(C, (mean(Xca) - mean(Xco)) / sqrt(vcase / ncase1 + vcont / ncont1))
				}## enf for j
			}
			else {
				for (j in 1 : nboot) {
					k = sample(length(X));
					Xnew = X[k]
						X1 = Xnew[Y == 1]; X2 = Xnew[Y == 0]
						C = c(C, calc_score_test(X1, X2))
				} ## end for j
			}## else if RVS
				cc = (sum(abs(C) >= abs(Tobs)) + 1) / (nboot + 1)
						return(cc)
}
*/