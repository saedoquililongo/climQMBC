function prob = getCDF(pdf,series,mu,sigma,skew,skewy)
%% getCDF:
%  This function evaluates each row of the series in the respective
%  cumulative distribution function assigned by the Kolmogorov-Smirnov (KS)
%  test in the getDist function of the climQMBC package. The available
%  distributions are:
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
% Input:
%   pdf = A vector with an ID for the resulting distribution from the KS
%         test. The ID is related to  the numeration of the distribution
%         listed in the description of this function.
%
%   series = An array of daily, monthly or annual data, without considering
%            leap days. Possible input dimensions:
%            2D array: [sub-periods, years]
%                      [sub-periods + days window, years]
%
%   mu = A vector with mean values of each sub-period.
%
%   sigma =  A vector with standard deviation values of each sub-period.
%
%   skew = A vector with skewness values of each sub-period.
%
%   skewy = A vector with skewness values of the logarithm of the series of
%           each sub-period.
%
% Output:
%   prob = A matrix with the non-exceedance probability for value row of
%          the input series.
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
% Get the number of rows and columns of the series.
rows = size(series,1);
cols = size(series,2);

% Compute the cumulative distribution function to the values of each row,
% based on the distribution assigned in the getDist function of the
% climQMBC package.
prob = zeros(rows,cols);
for sp=1:rows
    if pdf(sp) == 1 % i) Normal distribution.
        prob(sp,:) = normcdf(series(sp,:),mu(sp),sigma(sp));
        
    elseif pdf(sp) == 2 % ii) Log-Normal distribution.
        sigmay = sqrt(log(1+(sigma(sp)/mu(sp))^2));
        muy = log(mu(sp))-(sigmay^2)/2;
        prob(sp,:) = logncdf(series(sp,:),muy,sigmay);
        
    elseif pdf(sp) == 3 % iii) Gamma 2 parameters distribution.
        A = (sigma(sp)^2)/mu(sp);
        B = (mu(sp)/sigma(sp))^2;

        prob(sp,:) = gamcdf(series(sp,:),B,A);
        
    elseif pdf(sp) == 4 % iv) Gamma 3 parameters distribution.
        Bet  = (2/skew(sp))^2;
        Alp  = sigma(sp)/sqrt(Bet);
        Gam  = mu(sp)-(Alp*Bet);
        prob(sp,:)  = gamcdf((series(sp,:)-Gam),Bet,Alp); 
        
    elseif pdf(sp) == 5 % v) Log-Gamma 3 parameters distribution.
        Bety = (2/skewy(sp))^2;
        sigmay = sqrt(log(1+(sigma(sp)/mu(sp))^2));
        Alpy = sigmay/sqrt(Bety);
        muy = log(mu(sp))-(sigmay^2)/2;
        Gamy = muy - (Alpy*Bety);
        Lnsortdata = log(series(sp,:));
        Lnsortdata(isinf(Lnsortdata))= log(0.01);
        Lnsortdata(imag(Lnsortdata)~=0)= log(0.01); 
        prob(sp,:)  = gamcdf((Lnsortdata-Gamy),Bety,Alpy);
        
    elseif pdf(sp) == 6 % vi) Gumbel distribution.
        Sn = pi/sqrt(6);
        yn = 0.5772;
        a = Sn/sigma(sp);
        u = mu(sp)-(yn/a);
        prob(sp,:) = exp(-exp(-a.*(series(sp,:)-u)));
        
    elseif pdf(sp) == 7 % vii) Exponential distribution.
        gamexp = mu(sp) - sigma(sp);
        prob(sp,:) = max(1-exp(-1/sigma(sp)*(series(sp,:)-gamexp)),0);
    end
end

% Tail probabilities are set to (0 +  threshold) and (1 - threshold) to
% avoid numerical errors.
th = 10^-3;
prob(prob>1-th) = 1-th;
prob(prob<th) = th;
end