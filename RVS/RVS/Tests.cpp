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
	double q;
	int n;
	int ncontrol;
	double xbar;
	double xvar;
	double controlbar;
	double controlvar;

	for (size_t j = 0; j < snps.size(); j++) {
		p = 0;
		n = 0;
		ncontrol = 0;
		xbar = 0;
		xvar = 0;
		controlbar = 0;
		controlvar = 0;

		snp = snps[j];

		if (rvs) {
			for (size_t i = 0; i < snp.EG.size(); i++) {
				if (snp.EG[i] != NULL) {
					n++;
					p += IDmap[i];
					xbar += snp.EG[i];

					if (IDmap[i] == 0) {
						controlbar += snp.EG[i];
						ncontrol++;
					}
				}
			}

			xbar /= n;
			controlbar /= ncontrol;

			p /= n;
			q = 1 - p;

			for (size_t i = 0; i < snp.EG.size(); i++) {
				if (snp.EG[i] != NULL) {
					xvar += pow(snp.EG[i] - xbar, 2);

					if (IDmap[i] == 0)
						controlvar += pow(snp.EG[i] - controlbar, 2);
				}
			}

			xvar /= n - 1;
			controlvar /= ncontrol - 1;
			v = q * calcRobustVar(snp.p) + p * controlvar;
			v = xvar / v;
		}
		else
			v = 1;

		pvals.push_back(chiSquareOneDOF(scoreTest(IDmap, snp.EG) * v));
	}

	return pvals;
}
