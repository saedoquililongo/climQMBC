function [pdf,ks_fail] = getDist(series,allow_negatives, mu,sigma,skew,skewy)
%% getDist:
%  This function assigns an independent probability distribution function
%  to each row of the input series by comparing the empirical probability
%  distribution function with seven distributions based on the
%  Kolmogorov-Smirnov (KS) test. The available distributions are:
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
%   For allow_negatives=0, only 2), 3) and 5) are considered (1, 4, 6, and
%   7 are discarded). For series with negative values, only 1), 3), 4), 6),
%   and 7) are considered (2, 3 and 5 are discarded).
%
% Input:
%   series = An array of daily, monthly or annual data, without considering
%            leap days. Possible input dimensions:
%            2D array: [sub-periods, years]
%                      [sub-periods + days window, years]
%
%   allow_negatives = A flag that identifies if data allows negative values
%                     and also to replace no-rain values with random small 
%                     values (Chadwick et al., 2023) to avoid numerical
%                     problems with the probability distribution functions.
%                     allow_negatives = 1 or True: Allow negatives
%                     allow_negatives = 0 or False: Do not allow negative
%
%   mu = A vector with mean values of each sub-period.
%
%   sigma = A vector with standard deviation values of each sub-period.
%
%   skew = A vector with skewness values of each sub-period.
%
%   skewy = A vector with skewness values of the logarithm of the series
%           of each sub-period.
%
% Output:
%   pdf = A vector with an ID for the resulting distribution from the KS
%         test. The ID is related to  the numeration of the distribution
%         listed in the description of this function. This ID is used in
%         the getCDF and getCDFinv  functions of the climQMBC package.
%
% Written by Sebastian Aedo Quililongo (1*)
%            Cristian Chadwick         (2)
%            Fernando Gonzalez-Leiva   (3)
%            Jorge Gironas             (3, 4)
%           
%  (1) Stockholm Environment Institute, Latin America Centre, Bogota,
%      Colombia
%  (2) Faculty of Engineering and Sciences, Universidad Adolfo Ibanez,
%      Santiago, Chile
%  (3) Department of Hydraulics and Environmental Engineering, Pontificia 
%      Universidad Catolica de Chile, Santiago, Chile
%  (4) Centro de Cambio Global UC, Pontificia Universidad Catolica de Chile,
%      Santiago, Chile
%
% *Maintainer contact: sebastian.aedo.q@gmail.com
% Revision: 1, updated Apr 2024
%

%%
% Get the number of rows of the input series.
nrows = size(series,1);

% Initialize column vectors for the statistics needed for the available
% probability distribution functions.
pdf = zeros(nrows,1);

% Initialize a counter of failures of the KS test
ks_fail = 0;

% Perform the Kolmogorov-Smirnov test for sub-period or row.
for sp=1:nrows
    series_sub = series(sp,:);

    % Get the number of years with valid data (for allow_negatives=0, the 
    % series might have nan values to ignore them)
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
    
    % c) allow_negatives=0, set KS=1 to distributions that allow negative
    %    values (this will discard those distributions).
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