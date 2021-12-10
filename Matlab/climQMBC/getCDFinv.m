function xhat = getCDFinv(PDF,Taot,mu,sigma,skew,skewy)
%% getCDFinv:
%   This function evaluates the probability, Taot, in the respective
%   inverse cumulative distribution function assigned by the
%   Kolmogorov-Smirnov (KS) test in the getDist function of the climQMBC
%   package. The available distributions are:
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
%   Note that each row is associated to an independent distribution
%   function.
%
% Description:
%   1) Get the number of rows and years of the series.
%
%   2) Compute the inverse cumulative distribution function to the values
%      of each row, based on the distribution assigned in the getDist
%      function of the climQMBC package.
%
% Input:
%   PDF = A column vector with an ID for the resulting distribution from
%         the KS test [12,1]. The ID is related to the order of the 
%         distribution listed in the description of this function.
%
%   Taot = A column vector with the non-exceedance probability for each row
%          of the input series. [12,1] if the series consider monthly data
%          and [1,1] if the series consider annual data.
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
%   xhat = A column vector with the values obtained when the inverse
%          cumulative distribution function is applied. [12,1] if the
%          series consider monthly data and [1,1] if the series consider
%          annual data.
%

% Written by Sebastian Aedo Quililongo (1)
%            Cristian Chadwick         (2)
%            Fernando Gonzalez-Leiva   (1)
%            Jorge Gironas             (1)
%            
%   (1) Pontificia Universidad Catolica de Chile, Santiago, Chile
%       Department of Environmental and Hydraulic Engineering
%   (2) Universidad Adolfo Ibanez, Santiago, Chile
%       Faculty of Engineering and Sciences
% Maintainer contact: slaedo@uc.cl
% Revision: 0, updated Dec 2021


%%
% 1) Get the number of rows and years of the series.
n_m = size(Taot,1);
n_y = size(Taot,2);

% 2) Compute the inverse cumulative distribution function to the values of
%    each row, based on the distribution assigned in the getDist function
%    of the climQMBC package.
xhat = zeros(n_m,n_y);
for m=1:n_m
    if PDF(m)==1 % i) Normal distribution.
        xhat(m,:)  = norminv(Taot(m,:),mu(m),sigma(m));
        
    elseif PDF(m)==2 % ii) Log-Normal distribution.
        sigmay = sqrt(log(1+(sigma(m)/mu(m))^2));
        muy    = log(mu(m))-(sigmay^2)/2;
        xhat(m,:)  = logninv(Taot(m,:),muy,sigmay);
        
    elseif PDF(m)==3 % iii) Gamma 2 parameters distribution.
            A  = (sigma(m)^2)/mu(m);
            B  = (mu(m)/sigma(m))^2;
            
        xhat(m,:)= gaminv(Taot(m,:),B,A);
        
    elseif PDF(m)==4 % iv) Gamma 3 parameters distribution.
        Bet = (2/skew(m))^2; 
        Alp = sigma(m)/sqrt(Bet);
        Gam = mu(m)-(Alp*Bet);
        xhat(m,:)= (gaminv(Taot(m,:),Bet,Alp)+ Gam);
        
    elseif PDF(m)==5 % v) Log-Gamma 3 parameters distribution.
        BetyAster  = (2/skewy(m))^2;
        sigmay= sqrt(log(1+(sigma(m)/mu(m))^2));
        Alpy  = sigmay/sqrt(BetyAster);
        muy   = log(mu(m))-(sigmay^2)/2;
        Gamy  = muy -(Alpy*BetyAster);
        xhat(m,:)= exp(gaminv(Taot(m,:),BetyAster,Alpy)+Gamy);
        
    elseif PDF(m)==6 % vi) Gumbel distribution.
        Sn     = pi/sqrt(6);
        yn = 0.5772;
        a = Sn/sigma(m);
        u = mu(m)-(yn/a);
        xhat(m,:)= (u-log(-log(Taot(m,:)))./a);
        
    elseif PDF(m)==7 % vii) Exponential distribution.
        gamexp = mu(m) -sigma(m);
        xhat(m,:)= (gamexp-(sigma(m)*log(1-Taot(m,:))));           
    end     
end
end