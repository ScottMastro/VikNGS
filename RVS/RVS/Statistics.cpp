#include "stdafx.h"
#include "RVS.h"

/*
Uses EM algorithm to estimate the genotype frequencies in the sample

@param snp A SNP with genotype likelihoods.
@return A vector with three doubles to stands for probability of 0, 1 or 2 minor alleles.
*/
std::vector<double> calcEM(SNP &snp) {
	double p = 0.15;
	double q = 0.15;
	double qn = 1;
	double pn = 0;
	double dn = 0;
	double d = 0;

	int glCounter = 0;

	double Ep;
	double Eq;
	double pD;

	int k = 0;

	while (pow((pn - p), 2) + pow((qn - q), 2) > 0.000001) {

		d = 1 - p - q;

		glCounter = 0;
		Ep = 0;
		Eq = 0;

		for (size_t i = 0; i < snp.gl.size(); i++) {
			if (snp.L00(i) == NULL)
				continue;

			glCounter++;

			pD = 1 / (p * snp.L00(i) + q * snp.L01(i) + d * snp.L11(i));
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

	std::vector<double> freq;
	freq.push_back(p);
	freq.push_back(q);
	freq.push_back(1 - p - q);

	return freq;
}


/*
Function to calculate the conditional expected genotype probability E( P(G_ij | D_ij) ) 
given the genotype likelihoods P(D_ij | G_ij = g) and frequencies.

E( P(G_ij | D_ij) ) = sum from g=0 to 2 P(G_ij = g | D_ij),
where P(G_ij = g | D_ij) = P(D_ij | G_ij = g) * P(G_ij = g)/P(D_ij).

@param snp SNP with genotype likelihoods P(D_ij | G=AA, Aa or aa)} for one locus.
@return A vector containing conditional expectation probability E( P(G_ij | D_ij) ).
*/
std::vector<double> calcEG(SNP &snp) {

	double m0;
	double m1;
	double m2;
	double m;

	std::vector<double> EG;

	for (size_t i = 0; i < snp.gl.size(); i++) {

		if (snp.L00(i) != NULL) {
			m0 = snp.L00(i) * snp.p[0];
			m1 = snp.L01(i) * snp.p[1];
			m2 = snp.L11(i) * snp.p[2];
			m = 1 / (m0 + m1 + m2);

			EG.push_back(m1*m + 2 * m2*m);
		}
		else
			EG.push_back(NULL);

	}

	return EG;
}
