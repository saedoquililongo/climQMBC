function [pdf,ks_fail] = getDist(series,allow_negatives, mu,sigma,skew,skewy)
%% getDist:
%   This function assigns an independent probability distribution function
%   to each row of the input series by comparing the empirical probability
%   distribution function with seven distributions based on the
%   Kolmogorov-Smirnov (KS) test. If the series consider monthly data, it
%   will have 12 rows and each row will represent a month. For annual
%   data the series will have only one row. Only strictly positive
%   distributions are considered for precipitation and strictly positive
%   distributions are discarded if the series has negative values. The
%   available distributions are:
%          1) Normal distribution
%          2) Log-Normal distribution
%          3) Gamma 2 parameters distribution
%          4) Gamma 3 parameters distribution
%             (Pearson 3 parameters distribution)
%          5) Log-Gamma 3 parameters distribution
%             (Log-Pearson 3 parameters distribution)
%          6) Gumbel distribution
%          7) Exponential distribution
%
%   For precipitation, only 2), 3) and 5) are considered (1, 4, 6, and 7
%   are discarded). For series with negative values, only 1), 3), 4), 6),
%   and 7) are considered (2 and 5 are discarded).
%
% Description:
%   1) Get the number of years to compute the empirical distribution in
%      step 3) and get the number of rows of the input series.
%
%   2) Initialize column vectors for the statistics needed for the
%      available probability distribution functions.
%
%   3) Perform the Kolmogorov-Smirnov test for each row.
%      a) Get empirical distribution.
%      b) Compare distributions.
%      c) If variable is precipitation, discard distribution that allows 
%         negative values.
%      d) Compare KS values.
%
% Input:
%   series = A matrix of monthly or annual data (temperature or 
%            precipitation). If the series consider monthly data, it will
%            have 12 rows and each row will represent a month. For annual
%            data the series will have only one row.
%
%   mu = A column vector of mean values of the series. [12,1] if the series
%        consider monthly data and [1,1] if the series consider annual
%        data.
%
%   sigma = A column vector of standard deviation of the series. [12,1] if 
%           the series consider monthly data and [1,1] if the series
%           consider annual data.
%
%   skew = A column vector of skewness of the series. [12,1] if the series
%          consider monthly data and [1,1] if the series consider annual
%          data.
%
%   skewy = A column vector of skewness of the logarithm of the series. 
%           [12,1] if the series consider monthly data and [1,1] if the
%           series consider annual data.
%
%   var = A flag that identifies if data are temperature or precipitation.
%         This flag tells the getDist function if it has to discard
%         distribution functions that allow negative numbers.
%         Temperature:   var = 0
%         Precipitation: var = 1
%
% Output:
%   PDF = A column vector with an ID for the resulting distribution from
%         the KS test [12,1]. The ID is related to the order of the 
%         distribution listed in the description of this function. This ID
%         is used in the getCDF and getCDFinv functions of the climQMBC
%         package.
%

% Written by Sebastian Aedo Quililongo (1*)
%            Cristian Chadwick         (2)
%            Fernando Gonzalez-Leiva   (3)
%            Jorge Gironas             (3)
%            
%   (1) Centro de CAmbio Global UC, Pontificia Universidad Catolica de 
%       Chile, Santiago, Chile
%   (2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez,
%       Santiago, Chile
%   (3) Department of Hydraulics and Environmental Engineering, Pontificia
%       Universidad Catolica de Chile, Santiago, Chile
%       
%   *Maintainer contact: slaedo@uc.cl
% Revision: 0, updated Dec 2021

%%
% 1) Get the number of years to compute the empirical distribution in step 
%    3) and get the number of rows of the input series.
nrows = size(series,1);

% 2) Initialize column vectors for the statistics needed for the available
%    probability distribution functions.
pdf = zeros(nrows,1);

ks_fail = 0;
% 3) Perform the Kolmogorov-Smirnov test for each row.
for sp=1:nrows
    series_sub = series(sp,:);

    series_sub = series_sub(~isnan(series_sub));
    y_series = length(series_sub);

    % a) Get empirical distribution.
    sortdata = sort(series_sub)'; 
    probEmp = (1/(y_series+1):1/(y_series+1):y_series/(y_series+1))';

    % b) Compare distributions.
    % i) Normal distribution.
    normal = normcdf(sortdata,mu(sp),sigma(sp));
    KSnormal = max(abs(normal-probEmp));

    % ii) Log Normal distribution.
    if any(series(sp,:) < 0)
        KSlognormal = 1;
    else
        sigmay= sqrt(log(1+(sigma(sp)/mu(sp))^2));
        muy   = log(mu(sp))-(sigmay^2)/2;
        lognormal  = logncdf(sortdata,muy,sigmay);
        KSlognormal = max(abs(lognormal-probEmp));
    end
    
    % iii) Gamma 2 parameters distribution.
    if any(series(sp,:) < 0)
        KSgammaII = 1;
    else
        A   = (sigma(sp)^2)/mu(sp); 
        B   = (mu(sp)/sigma(sp))^2;
        GammaII  = gamcdf(sortdata,B,A);
        KSgammaII = max(abs(GammaII-probEmp));
    end

    % iv) Gamma 3 parameters distribution.
    %     (Pearson 3 parameters distribution)
    Bet  = (2/skew(sp))^2; 
    Alp  = sigma(sp)/sqrt(Bet);
    Gam  = mu(sp)-(Alp*Bet);
    GammaIII  = gamcdf((sortdata-Gam),Bet,Alp);
    KSgammaIII = max(abs(GammaIII-probEmp));

    % v) Log-Gamma 3 parameters distribution.
    %    (Log-Pearson 3 parameters distribution)
    if any(series(sp,:) < 0)
        KSLpIII = 1;
    else
        Bety = (2/skewy(sp))^2;
        Alpy = sigmay/sqrt(Bety);
        Gamy = muy-(Alpy*Bety);
        Lnsortdata = log(sortdata);
        Lnsortdata(isinf(Lnsortdata)) = log(0.01); 
        LpIII  = gamcdf((Lnsortdata-Gamy),Bety,Alpy);
        KSLpIII = max(abs(LpIII-probEmp));
    end

    % vi) Gumbel distribution.
    Sn    = pi/sqrt(6);
    yn    = 0.5772;
    a = Sn/sigma(sp);
    u = mu(sp)-(yn/a);
    gumbel = exp(-exp(-a.*(sortdata-u)));
    KSgumbel = max(abs(gumbel-probEmp));

    % vii) Exponential distribution.
    gamexp = mu(sp)-sigma(sp);
    exponential = max(1-exp(-1/sigma(sp)*(sortdata-gamexp)),0);
    KSexponential= max(abs(exponential-probEmp));
    
    % c) If variable is precipitation, set KS=1 to distributions that allow
    %    negative values (this will discard those distributions).
    if allow_negatives==0
        KSnormal=1;
        KSgammaIII=1;
        KSgumbel=1;
        KSexponential=1;

        KSLpIII = 1;
    end
    
    % d) The distribution with lower KS value is considered for each month.
    KS_vals = [KSnormal,KSlognormal,KSgammaII,KSgammaIII, KSLpIII,KSgumbel,KSexponential];
    [min_KS,bestPDF]= min(KS_vals);

    ks_crit = 1.3581/sqrt(y_series);
    ks_fail = ks_fail + (sign(min_KS-ks_crit)+1)/2;

    pdf(sp)= bestPDF;
end
end