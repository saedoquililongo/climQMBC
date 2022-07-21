function PDF = getDist(series,mu,sigma,skew,skewy,var)
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
y_series = size(series,2);
n = size(series,1);

% 2) Initialize column vectors for the statistics needed for the available
%    probability distribution functions.
PDF = zeros(n,1);
sigmay = zeros(n,1);
muy    = zeros(n,1);
A      = zeros(n,1);
B      = zeros(n,1);
Alp    = zeros(n,1);
Bet    = zeros(n,1);
Gam    = zeros(n,1);
Alpy   = zeros(n,1);
Bety   = zeros(n,1);
Gamy   = zeros(n,1);
a      = zeros(n,1);
u      = zeros(n,1);

% 3) Perform the Kolmogorov-Smirnov test for each row.
for m=1:n
    
    % a) Get empirical distribution.
    sortdata = sort(series(m,:))'; 
    probEmp = (1/(y_series+1):1/(y_series+1):y_series/(y_series+1))';

    % b) Compare distributions.
    % i) Normal distribution.
    normal = normcdf(sortdata,mu(m),sigma(m));
    KSnormal = max(abs(normal-probEmp));

    % ii) Log Normal distribution.
    if any(series(m,:) < 0)
        KSlognormal = 1;
    else
        sigmay(m)= sqrt(log(1+(sigma(m)/mu(m))^2));
        muy(m)   = log(mu(m))-(sigmay(m)^2)/2;
        lognormal  = logncdf(sortdata,muy(m),sigmay(m));
        KSlognormal = max(abs(lognormal-probEmp));
    end
    
    % iii) Gamma 2 parameters distribution.
    if any(series(m,:) < 0)
        KSgammaII = 1;
    else
        A(m)   = (sigma(m)^2)/mu(m); 
        B(m)   = (mu(m)/sigma(m))^2;
        GammaII  = gamcdf(sortdata,B(m),A(m));
        KSgammaII = max(abs(GammaII-probEmp));
    end

    % iv) Gamma 3 parameters distribution.
    %     (Pearson 3 parameters distribution)
    Bet(m)  = (2/skew(m))^2; 
    Alp(m)  = sigma(m)/sqrt(Bet(m));
    Gam(m)  = mu(m)-(Alp(m)*Bet(m));
    GammaIII  = gamcdf((sortdata-Gam(m)),Bet(m),Alp(m));
    KSgammaIII = max(abs(GammaIII-probEmp));

    % v) Log-Gamma 3 parameters distribution.
    %    (Log-Pearson 3 parameters distribution)
    if any(series(m,:) < 0)
        KSLpIII = 1;
    else
        Bety(m) = (2/skewy(m))^2;
        Alpy(m) = sigmay(m)/sqrt(Bety(m));
        Gamy(m) = muy(m)-(Alpy(m)*Bety(m));
        Lnsortdata = log(sortdata);
        Lnsortdata(isinf(Lnsortdata)) = log(0.01); 
        LpIII  = gamcdf((Lnsortdata-Gamy(m)),Bety(m),Alpy(m));
        KSLpIII = max(abs(LpIII-probEmp));
    end

    % vi) Gumbel distribution.
    Sn    = pi/sqrt(6);
    yn    = 0.5772;
    a(m) = Sn/sigma(m);
    u(m) = mu(m)-(yn/a(m));
    gumbel = exp(-exp(-a(m).*(sortdata-u(m))));
    KSgumbel = max(abs(gumbel-probEmp));

    % vii) Exponential distribution.
    gamexp = mu(m)-sigma(m);
    exponential = max(1-exp(-1/sigma(m)*(sortdata-gamexp)),0);
    KSexponential= max(abs(exponential-probEmp));
    
    % c) If variable is precipitation, set KS=1 to distributions that allow
    %    negative values (this will discard those distributions).
    if var==1
        KSnormal=1;
        KSgammaIII=1;
        KSgumbel=1;
        KSexponential=1;
    end
    
    % d) The distribution with lower KS value is considered for each month.
    [~,bestPDF]= min([KSnormal,KSlognormal,KSgammaII,KSgammaIII,...
        KSLpIII,KSgumbel,KSexponential]);

    PDF(m)= bestPDF;
end
end