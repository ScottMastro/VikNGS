#include "stdafx.h"
#include "RVS.h"

#include <iostream>  
#include <math.h> 

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

std::vector<double> LogisticRegression(std::vector<bool> &y, std::vector<double> &x) {

	double w = 0;
	double w0 = 0;
	double learningRate = 0.02;

	double wlast = 100;
	double w0last = 100;

	double wgrad;
	double w0grad;

	double tmp = 0;
	int iteration = 0;

	while (true) {
		iteration++;
		wgrad = 0;
		w0grad = 0;
		for (size_t i = 0; i < y.size(); i++) {
			if (x[i] != NULL) {
				tmp = y[i] - (1.0 / (1.0 + exp(-(w * x[i] + w0))));
				wgrad += x[i] * tmp;
				w0grad += tmp;
			}
		}

		w += learningRate * wgrad;
		w0 += learningRate * w0grad;

		if (iteration > 1000 || (abs(w - wlast) < 0.0000001 && abs(w0 - w0last) < 0.0000001)) {
			std::vector<double> out;
			out.push_back(w0);
			out.push_back(w);
			std::cout << w0;
			std::cout << '\n';

			return out;
		}

		wlast = w;
		w0last = w0;
	}
}

double LogisticRegressionInterceptOnly(std::vector<bool> &y, std::vector<double> &x) {

	double w0 = 0;
	double learningRate = 0.02;

	double w0last = 100;
	double w0grad;

	int iteration = 0;

	while (true) {
		iteration++;
		w0grad = 0;
		for (size_t i = 0; i < y.size(); i++)
			if (x[i] != NULL)
				w0grad += y[i] - (1.0 / (1.0 + exp(-w0)));

		w0 += learningRate * w0grad;

		if (iteration > 1000 || abs(w0 - w0last) < 0.0000001)
			return w0;

		w0last = w0;

		std::cout << w0;
		std::cout << '\n';
	}
}