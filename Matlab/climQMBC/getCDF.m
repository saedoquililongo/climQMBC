function Taot = getCDF(PDF,series,mu,sigma,skew,skewy)
%% getCDF:
%   This function evaluates each row of the series in the respective
%   cumulative distribution function assigned by the Kolmogorov-Smirnov
%   (KS) test in the getDist function of the climQMBC package. The
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
%   Note that each row is associatied to an independent distribution
%   function.
%
% Description:
%   1) Get the number of rows and years of the series.
%   
%   2) Compute the cumulative distribution function to the values of each 
%      row, based on the distribution assigned in the getDist function of
%      the climQMBC package.
%
%   3) Tail probabilities are set to (0 + threshold) and (1 - threshold) to
%      avoid numerical errors. Threshold is set in the last lines of this
%      function. Default is th = 10^-3.
%
% Input:
%   PDF = A column vector with an ID for the resulting distribution from
%         the KS test [12,1]. The ID is related to the order of the 
%         distribution listed in the description of this function.
%
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
% Output:
%   Taot = A column vector with the non-exceedance probability for each row
%          of the input series. [12,1] if the series consider monthly data
%          and [1,1] if the series consider annual data.
%

% Written by Sebastian Aedo Quililongo (1)
%            Cristian Chadwick         (2)
%            Fernando Gonzalez-Leiva   (1)
%            Jorge Gironas             (1)
%            
%   (1) Pontificia Universidad Catolica de Chile, Santiago, Chile
%   (2) Universidad Adolfo Ibanez, Santiago, Chile
%       Faculty of Engineering and Sciences
%       Faculty of Forest Science
% Maintainer contact: slaedo@uc.cl
% Revision: 0, updated Dec 2021


%%
% 1) Get the number of rows and years of the series.
n_m = size(series,1);
n_y = size(series,2);

% 2) Compute the cumulative distribution function to the values of each 
%    row, based on the distribution assigned in the getDist function of the
%    climQMBC package.
Taot = zeros(n_m,n_y);
for m=1:n_m
    if PDF(m) == 1 % i) Normal distribution.
        Taot(m,:) = normcdf(series(m,:),mu(m),sigma(m));
        
    elseif PDF(m) == 2 % ii) Log-Normal distribution.
        sigmay = sqrt(log(1+(sigma(m)/mu(m))^2));
        muy = log(mu(m))-(sigmay^2)/2;
        Taot(m,:) = logncdf(series(m,:),muy,sigmay);
        
    elseif PDF(m) == 3 % iii) Gamma 2 parameters distribution.
        A = (sigma(m)^2)/mu(m);
        B = (mu(m)/sigma(m))^2;

        Taot(m,:) = gamcdf(series(m,:),B,A);
        
    elseif PDF(m) == 4 % iv) Gamma 3 parameters distribution.
        Bet  = (2/skew(m))^2;
        Alp  = sigma(m)/sqrt(Bet);
        Gam  = mu(m)-(Alp*Bet);
        Taot(m,:)  = gamcdf((series(m,:)-Gam),Bet,Alp); 
        
    elseif PDF(m) == 5 % v) Log-Gamma 3 parameters distribution.
        Bety = (2/skewy(m))^2;
        sigmay = sqrt(log(1+(sigma(m)/mu(m))^2));
        Alpy = sigmay/sqrt(Bety);
        muy = log(mu(m))-(sigmay^2)/2;
        Gamy = muy - (Alpy*Bety);
        Lnsortdata = log(series(m,:));
        Lnsortdata(isinf(Lnsortdata))= log(0.01);
        Lnsortdata(imag(Lnsortdata)~=0)= log(0.01); 
        Taot(m,:)  = gamcdf((Lnsortdata-Gamy),Bety,Alpy);
        
    elseif PDF(m) == 6 % vi) Gumbel distribution.
        Sn = pi/sqrt(6);
        yn = 0.5772;
        a = Sn/sigma(m);
        u = mu(m)-(yn/a);
        Taot(m,:) = exp(-exp(-a.*(series(m,:)-u)));
        
    elseif PDF(m) == 7 % vii) Exponential distribution.
        gamexp = mu(m) -sigma(m);
        Taot(m,:) = max(1-exp(-1/sigma(m)*(series(m,:)-gamexp)),0);
    end
end

% 3) Tail probabilities are set to (0 +  threshold) and (1 - threshold) to
%    avoid numerical errors.
th = 10^-3;
Taot(Taot>1-th) = 1-th;
Taot(Taot<th) = th;
end