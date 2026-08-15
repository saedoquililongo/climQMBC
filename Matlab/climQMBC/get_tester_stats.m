% This script contains computation of the tester statistics, including the 
% NSE, KGE and ratios (or difference) in the mean, std. dev., skewness and
% several percentiles.
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
% Revision: 0, updated Aug 2026

function [li_stats] = get_tester_stats(bc_series,obj_series,mult_change)
bc_mean = mean(bc_series);
obj_mean = mean(obj_series);

a =  sum((bc_series-bc_mean).*(obj_series-obj_mean));
b = sqrt(sum((bc_series-bc_mean).^2));
c = sqrt(sum((obj_series-obj_mean).^2));

r2 = (a/(b*c))^2;

    if mult_change == 1
        mu_ratio = mean(bc_series)/mean(obj_series);
        sigma_ratio = std(bc_series)/std(obj_series);
        sk_ratio = skewness(bc_series)/skewness(obj_series);
        p05_ratio = prctile(bc_series,5)/prctile(obj_series,5);
        p10_ratio = prctile(bc_series,10)/prctile(obj_series,10);
        p25_ratio = prctile(bc_series,25)/prctile(obj_series,25);
        p50_ratio = prctile(bc_series,50)/prctile(obj_series,50);
        p75_ratio = prctile(bc_series,75)/prctile(obj_series,75);
        p90_ratio = prctile(bc_series,90)/prctile(obj_series,90);
        p95_ratio = prctile(bc_series,95)/prctile(obj_series,95);

        kge = 1-sqrt((1-r2)^2 + (1-mu_ratio)^2 + (1-sigma_ratio)^2);
        nse = 1-sum((bc_series-obj_series).^2)/sum((mean(obj_series)-obj_series).^2);
    else
        mu_ratio = mean(bc_series)-mean(obj_series);
        sigma_ratio = std(bc_series)-std(obj_series);
        sk_ratio = skewness(bc_series)-skewness(obj_series);
        p05_ratio = prctile(bc_series,5)-prctile(obj_series,5);
        p10_ratio = prctile(bc_series,10)-prctile(obj_series,10);
        p25_ratio = prctile(bc_series,25)-prctile(obj_series,25);
        p50_ratio = prctile(bc_series,50)-prctile(obj_series,50);
        p75_ratio = prctile(bc_series,75)-prctile(obj_series,75);
        p90_ratio = prctile(bc_series,90)-prctile(obj_series,90);
        p95_ratio = prctile(bc_series,95)-prctile(obj_series,95);

        kge = sqrt((1-r2)^2 + (1-mu_ratio)^2 + (1-sigma_ratio)^2);
        nse = sum((bc_series-obj_series).^2)/sum((mean(obj_series)-obj_series).^2);
    end
    
    li_stats = [nse,kge,mu_ratio,sigma_ratio,sk_ratio,p05_ratio,p10_ratio,p25_ratio,p50_ratio,p75_ratio,p90_ratio,p95_ratio];
end