function [eta_cvtd,beta_cvtd] = frag_conversion(point1,point2)
%{
This is the code for the fragility conversion method by Pang and Wang (202X):

Reference: 
Pang Y., and Wang X. (202X) “Cloud-IDA-MSA conversion of fragility curves
for efficient and high-fidelity resilience assessment”. Journal of
Structural Engineering (ASCE), under review (STENG-9556R1)

Input:
   point1 - 1-by-2 vector of the first converted fragility point (pf1, IM1)
   point2 - 1-by-2 vector of the second converted fragility point (pf2, IM2)

Output:
 eta_cvtd - converted fragility median (scaler)
beta_cvtd - converted fragility dispersion (scaler)
%}

IM1    = point1(2); IM2    = point2(2);
Pf_IM1 = point1(1); Pf_IM2 = point2(1);
eta_cvtd = exp((log(IM1)*norminv(Pf_IM2,0,1)-log(IM2)*norminv(Pf_IM1,0,1))...
    /(norminv(Pf_IM2,0,1)-norminv(Pf_IM1,0,1)));
beta_cvtd = (log(IM2)-log(IM1))...
    /(norminv(Pf_IM2,0,1)-norminv(Pf_IM1,0,1));

end

