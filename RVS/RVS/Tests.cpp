#include "stdafx.h"
#include "RVS.h"

#include <iostream>  
#include <vector>

void RVSasy(std::vector<SNP> &snps, std::vector<bool> &IDmap, bool rsv) {
	LogisticRegression(IDmap, snps[0].EG);
	LogisticRegressionInterceptOnly(IDmap, snps[0].EG);


}




/*
RVS using asymptotic distribution for score test statistic
#' 
#' This functions actually includes two association tests which differ only at the estimation of the variance of score statistic by option 'method'.
#' It uses regular formula for \eqn{Var_{case}(E(G_{ij}|D_{ij}))} when method = 'Regular' (also called likelihood method in the paper),
#' \eqn{Var_{case}(G_{ij})}  is used for \eqn{Var_{case}(E(G_{ij}|D_{ij}))} when method = 'RVS'. It uses the population frequency P
#' to estimate variance for case if method=RVS, otherwise use the expected genotype probability for case directly.
#'
#' Note, both function scaled the variance in equation (1) in Appendix A by dividing N_case*Ncontrol/N
#' also note the test statistics in anova gives us Rao=s^2/var(s^2), where RVS=s^2/robust var(s^2)
#' so RVS=Rao *var(s^2)/robust var(s^2), s=sum(y_j-mean(y))E(G_ij|D_ij), therefore var(s)//var(E(Gij|Dij)), var(X) in code.
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/22570057}
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/24733292}
#' @param Y - a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X - a vector of genotype for a snp, first ncase and then ncont
#' @param P (vector with length 3  (double)) - estimated genotype frequency for a variant P(G=0), P(G=1) and P(G=2)
#' @param RVS, logic variable to indicate the estimation of variance of the score statistic. 
#' @return   p-values for the snp
#' @export
RVS_asy = function(Y, X, P, RVS = 'TRUE') {

## run anlysis for non-NA X only

	Y = Y[!is.na(X)]
		X = X[!is.na(X)]
		a = (glm(Y~1, family = 'binomial'))
		b = glm(Y~X, family = 'binomial')
		p = length(X[Y == 1]) / length(X)
		q = 1 - p
#  if(RVS %in% c('TRUE','True','true','T','t')) {
		v = q*as.numeric(calc_robust_Var(P)) + p*var(X[Y == 0],na.rm = TRUE)
# }else{
		#   v = q*var(X[Y == 1]) + p*var(X[Y == 0])  ## regular estimation for variance on case at var(X[Y == 1]) from test_RVS_asy
# }
		x = anova(a, b, test = 'Rao')
		if (RVS %in% c('TRUE', 'True', 'true', 'T', 't')) {
			x_rvs = x$Rao[2] * var(X) / v  ## var(S_inRao) = #case*#cont / #sample*var(X), var(RVS) = #case*#cont / #sample*(q*var(Xcase) + p*var(Xcont))
		}
		else {
			x_rvs = x$Rao[2]
		}
		cc = 1 - pchisq(x_rvs, 1)
			# res<-list(coefficients = coef(summary(b))[2, 1],
				#       score = x$Rao[2],
				#      p_Rao = x[2, 6],
				#       RVS = x_rvs,
				#       p_rvs = cc)
			return (cc)
}
*/