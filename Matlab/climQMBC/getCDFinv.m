function xhat = getCDFinv(pdf,prob,mu,sigma,skew,skewy)
%% getCDFinv:
%  This function evaluates the probability, prob, in the respective inverse
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
%   Note that each row is associated to an independent distribution
%   function.
%
% Input:
%   pdf = A vector with an ID for the resulting distribution from the KS
%         test. The ID is related to  the numeration of the distribution
%         listed in the description of this function.
%
%   prob = A matrix with the non-exceedance probability for value row of
%          the input series.
%
%   mu = A vector with mean values of each sub-period.
%
%   sigma = A vector with standard deviation values of each sub-period.
%
%   skew = A vector with skewness values of each sub-period.
%
%   skewy = A vector with skewness values of the logarithm of the series of
%           each sub-period.
%
% Output:
%   xhat = A matrix with the values obtained when the inverse cumulative
%            distribution function is applied to prob.
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
rows = size(prob,1);
cols = size(prob,2);

% Compute the cumulative distribution function to the values of each row,
% based on the distribution assigned in the getDist function of the
% climQMBC package.
xhat = zeros(rows,cols);
for sp=1:rows
    if pdf(sp)==1 % i) Normal distribution.
        xhat(sp,:)  = norminv(prob(sp,:),mu(sp),sigma(sp));
        
    elseif pdf(sp)==2 % ii) Log-Normal distribution.
        sigmay = sqrt(log(1+(sigma(sp)/mu(sp))^2));
        muy    = log(mu(sp))-(sigmay^2)/2;
        xhat(sp,:)  = logninv(prob(sp,:),muy,sigmay);
        
    elseif pdf(sp)==3 % iii) Gamma 2 parameters distribution.
            A  = (sigma(sp)^2)/mu(sp);
            B  = (mu(sp)/sigma(sp))^2;
            
        xhat(sp,:)= gaminv(prob(sp,:),B,A);
        
    elseif pdf(sp)==4 % iv) Gamma 3 parameters distribution.
        Bet = (2/skew(sp))^2; 
        Alp = sigma(sp)/sqrt(Bet);
        Gam = mu(sp)-(Alp*Bet);
        xhat(sp,:)= (gaminv(prob(sp,:),Bet,Alp)+ Gam);
        
    elseif pdf(sp)==5 % v) Log-Gamma 3 parameters distribution.
        BetyAster  = (2/skewy(sp))^2;
        sigmay= sqrt(log(1+(sigma(sp)/mu(sp))^2));
        Alpy  = sigmay/sqrt(BetyAster);
        muy   = log(mu(sp))-(sigmay^2)/2;
        Gamy  = muy -(Alpy*BetyAster);
        xhat(sp,:)= exp(gaminv(prob(sp,:),BetyAster,Alpy)+Gamy);
        
    elseif pdf(sp)==6 % vi) Gumbel distribution.
        Sn     = pi/sqrt(6);
        yn = 0.5772;
        a = Sn/sigma(sp);
        u = mu(sp)-(yn/a);
        xhat(sp,:)= (u-log(-log(prob(sp,:)))./a);
        
    elseif pdf(sp)==7 % vii) Exponential distribution.
        gamexp = mu(sp) -sigma(sp);
        xhat(sp,:)= (gamexp-(sigma(sp)*log(1-prob(sp,:))));           
    end     
end
end