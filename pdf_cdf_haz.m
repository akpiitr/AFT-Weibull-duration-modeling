function [pdf, cdf, haz] = pdf_cdf_haz(t, CR, AADT)

cr = CR-1;
%lambda vals assuming AADT = 20000 below;
% lamda_vals = [0.229436265, 0.387688925, 0.305722332, 0.227036155, 0.176190761, 0.054554504, 0.038926112, 0.0839332873];
lambda_intercept = [2.73213, 1.547552, 1.625078, 1.722646, 1.936188, 2.148555, 2.68609, 1.597733];
lambda_AADT = [-0.000063, -0.00003, -0.000022, -0.000012, -0.00001, 0.000038, 0.000028, 0.000044];
p_vals = [1.902977948, 1.926303728, 1.964019228, 1.912218335, 1.814583564, 1.819500477, 1.908773711, 1.421540297];

lamda_vals =(1./exp(lambda_intercept+lambda_AADT*AADT));

lam = lamda_vals(cr);
p = p_vals(cr);

cdf = 1 - exp(-(lam*t)^p);

pdf = lam * p * exp(-(lam*t)^p) * (lam*t)^(p-1);

haz = lam*p*(lam*t)^(p-1);

end