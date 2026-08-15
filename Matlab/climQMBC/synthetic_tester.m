% This script contains the synthetic tester function to compare the performance 
% of different methods based on a synthetic example rather than the original
% data, to compare the theoretical performance.
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

function [tester_stats] = synthetic_tester(obs, mod, allow_negatives, mult_change, SDM_var, yn)
    if ~exist('yn','var')
      yn=10000;
    end

    % Get series statistics
    mu_obs_h = mean(obs);
    sigma_obs_h = std(obs);

    mu_mod_h = mean(mod(1:length(obs)));
    sigma_mod_h = std(mod(1:length(obs)));

    mu_mod_f = mean(mod(length(obs)+1:end));
    sigma_mod_f = std(mod(length(obs)+1:end));

    % Get bias and future values
    if mult_change == 1
        mu_bias = mu_mod_h/mu_obs_h - 1;
        sigma_bias = sigma_mod_h/sigma_obs_h - 1;

        mu_change = mu_mod_f/mu_mod_h - 1;
        sigma_change = sigma_mod_f/sigma_mod_h - 1;

        mu_obj_f = mu_obs_h*(1+mu_change);
        sigma_obj_f = sigma_obs_h*(1+sigma_change);

    else
        mu_bias = mu_mod_h-mu_obs_h;
        sigma_bias = sigma_mod_h-sigma_obs_h;

        mu_change = mu_mod_f-mu_mod_h;
        sigma_change = sigma_mod_f-sigma_mod_h;

        mu_obj_f = mu_obs_h + mu_change;
        sigma_obj_f = sigma_obs_h + sigma_change;
    end

    % Get synthetic series
    if allow_negatives == 0
        k = (mu_obs_h/sigma_obs_h)^2;
        th = (sigma_obs_h^2)/mu_obs_h;
        obs_h = gamrnd(k,th,yn,1);

        k = (mu_obj_f/sigma_obj_f)^2;
        th = (sigma_obj_f^2)/mu_obj_f;
        obs_f = gamrnd(k,th,yn,1);

    else
        obs_h = normrnd(mu_obs_h,sigma_obs_h,yn,1);
        obs_f = normrnd(mu_obj_f,sigma_obj_f,yn,1);
    end
    
    % Add bias considering an equiprobability (normal) transform
    mu_temp = mean(obs_h);
    sigma_temp = std(obs_h);

    if mult_change == 1
        mu_scaled = mu_temp*(1 + mu_bias);
        sigma_scaled = sigma_temp*(1 + sigma_bias);
    else
        mu_scaled = mu_temp + mu_bias;
        sigma_scaled = sigma_temp + sigma_bias;
    end

    mod_h = mu_scaled + (obs_h - mu_temp)*sigma_scaled/sigma_temp;
    mod_f = mu_scaled + (obs_f - mu_temp)*sigma_scaled/sigma_temp;


    % Concat historical and future period
    % Future period is concated twice so that only the second future is
    % analyzed, without moving windows
    obs = obs_h;
    mod = [mod_h;mod_f;mod_f];

    % Apply bias correction methods
    QM_h_series = QM(obs_h,mod_h,allow_negatives,'A');
    SDM_h_series = SDM(obs_h,mod_h,SDM_var,'A');
    QM_f_series = QM(obs,mod,allow_negatives,'A');
    QM_f_series = QM_f_series(2*yn+1:end);
    DQM_f_series = DQM(obs,mod,mult_change,allow_negatives,'A');
    DQM_f_series = DQM_f_series(2*yn+1:end);
    QDM_f_series = QDM(obs,mod,mult_change,allow_negatives,'A');
    QDM_f_series = QDM_f_series(2*yn+1:end);
    UQM_f_series = UQM(obs,mod,mult_change,allow_negatives,'A');
    UQM_f_series = UQM_f_series(2*yn+1:end);
    SDM_f_series = SDM(obs,mod,SDM_var,'A');
    SDM_f_series = SDM_f_series(2*yn+1:end);


    li_stats_qm_h = get_tester_stats(QM_h_series,obs_h,mult_change);
    li_stats_sdm_h = get_tester_stats(SDM_h_series,obs_h,mult_change);
    li_stats_qm_f = get_tester_stats(QM_f_series,obs_f,mult_change);
    li_stats_dqm_f = get_tester_stats(DQM_f_series,obs_f,mult_change);
    li_stats_qdm_f = get_tester_stats(QDM_f_series,obs_f,mult_change);
    li_stats_uqm_f = get_tester_stats(UQM_f_series,obs_f,mult_change);
    li_stats_sdm_f = get_tester_stats(SDM_f_series,obs_f,mult_change);

    tester_stats = [li_stats_qm_h;li_stats_sdm_h;li_stats_qm_f;li_stats_dqm_f;li_stats_qdm_f;li_stats_uqm_f;li_stats_sdm_f];
    rownames = {'QM_h','SDM_h','QM_f','DQM_f','QDM_f','UQM_f','SDM_f'};
    if mult_change == 1
        colnames = {'NSE','KGE','Mean','StdDev','Skew','P05','P10','P25','P50','P75','P90','P95'};
    else
        colnames = {'1 minus NSE','1 minus KGE','Mean','StdDev','Skew','P05','P10','P25','P50','P75','P90','P95'};
        colnames = matlab.lang.makeValidName(colnames);
    end
    tester_stats = array2table(tester_stats,'VariableNames',colnames,'RowNames',rownames);
end