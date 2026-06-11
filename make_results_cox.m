% make_results_cox.m
%
% Collects MCMC results for the log-Gaussian Cox process experiment
% and writes LaTeX table rows and log-likelihood trace figures.
%
% Methods included (Table 5):
%   mGrad, pMALA, Ellipt, pCN, pCNL, RHMC, GI-MALA
%
% Two configurations are tabulated (outer loop i = 1:2):
%   i = 1  ->  N = 64  (full resolution, ds = 1)  -- main Table 5
%   i = 2  ->  N = 32  (reduced resolution, ds = 2)
% RHMC and GI-MALA were run at N = 64 only; their rows appear in table 1.
%
% Output written to diagrams/:
%   GaussianCox_table1.txt   LaTeX rows for Table 5 (N=64)
%   GaussianCox_table2.txt   LaTeX rows for Table 5 (N=32)
%   GaussianCox_*_logL1.eps  Log-L traces for Figure 6
%
% Run from the repository root after run_all_experiments.m.

outdir = 'diagrams/';
addpath results/Cox_regression/;

datasetName = 'GaussianCox';
tableNames  = {'Table5_N64', 'Table5_N32'};
Repeats = 10;

for i = 1:2

    % Pre-allocate all ESS arrays with NaN for robustness
    [ESSminMarg,   ESSmedianMarg,   ESSmaxMarg,   ESSlogLMarg,   TrainTimeMarg,   deltaMarg,   accRateMarg  ] = deal(NaN(1,Repeats));
    [ESSminMALA,   ESSmedianMALA,   ESSmaxMALA,   ESSlogLMALA,   TrainTimeMALA,   deltaMALA,   accRateMALA  ] = deal(NaN(1,Repeats));
    [ESSminEllipt, ESSmedianEllipt, ESSmaxEllipt, ESSlogLEllipt, TrainTimeEllipt                            ] = deal(NaN(1,Repeats));
    [ESSminpCN,    ESSmedianpCN,    ESSmaxpCN,    ESSlogLpCN,    TrainTimepCN,    deltapCN,    accRatepCN   ] = deal(NaN(1,Repeats));
    [ESSminpCNL,   ESSmedianpCNL,   ESSmaxpCNL,   ESSlogLpCNL,   TrainTimepCNL,   deltapCNL,   accRatepCNL  ] = deal(NaN(1,Repeats));
    [ESSminRHMC,   ESSmedianRHMC,   ESSmaxRHMC,   ESSlogLRHMC,   TrainTimeRHMC                              ] = deal(NaN(1,Repeats));
    [ESSminGimala, ESSmedianGimala, ESSmaxGimala, ESSlogLGimala, TrainTimeGimala, deltaGimala, accRateGimala] = deal(NaN(1,Repeats));

    for rep = 1:Repeats

        % --- mGrad ---
        fname = ['logGaussianCoxGirolami_Marg_repeat' num2str(rep) '.mat'];
        if exist(fname,'file')
            load(fname);
            if iscell(summaryMarg) && length(summaryMarg) >= i
                ESSminMarg(rep)    = min(summaryMarg{i}.eff_F);
                ESSmedianMarg(rep) = median(summaryMarg{i}.eff_F);
                ESSmaxMarg(rep)    = max(summaryMarg{i}.eff_F);
                ESSlogLMarg(rep)   = summaryMarg{i}.eff_LogL;
                TrainTimeMarg(rep) = summaryMarg{i}.elapsed;
                deltaMarg(rep)     = summaryMarg{i}.delta;
                accRateMarg(rep)   = summaryMarg{i}.accRates.F;
            end
        end

        % --- pMALA (single config, i=1 only) ---
        if i == 1
            fname = ['logGaussianCoxGirolami_Marg_repeat' num2str(rep) '_Marg_pMALA.mat'];
            if exist(fname,'file')
                load(fname);
                ESSminMALA(rep)    = min(ess);
                ESSmedianMALA(rep) = median(ess);
                ESSmaxMALA(rep)    = max(ess);
                ESSlogLMALA(rep)   = eff_LogL;
                TrainTimeMALA(rep) = elapsedTime;
                deltaMALA(rep)     = delta;
                accRateMALA(rep)   = accRates.F;
            end
        end

        % --- Ellipt ---
        fname = ['logGaussianCoxGirolami_Ellipt_repeat' num2str(rep) '.mat'];
        if exist(fname,'file')
            load(fname);
            if iscell(summaryEllipt) && length(summaryEllipt) >= i
                ESSminEllipt(rep)    = min(summaryEllipt{i}.eff_F);
                ESSmedianEllipt(rep) = median(summaryEllipt{i}.eff_F);
                ESSmaxEllipt(rep)    = max(summaryEllipt{i}.eff_F);
                ESSlogLEllipt(rep)   = summaryEllipt{i}.eff_LogL;
                TrainTimeEllipt(rep) = summaryEllipt{i}.elapsed;
            end
        end

        % --- pCN ---
        fname = ['logGaussianCoxGirolami_pCN_repeat' num2str(rep) '.mat'];
        if exist(fname,'file')
            load(fname);
            if iscell(summarypCN) && length(summarypCN) >= i
                ESSminpCN(rep)    = min(summarypCN{i}.eff_F);
                ESSmedianpCN(rep) = median(summarypCN{i}.eff_F);
                ESSmaxpCN(rep)    = max(summarypCN{i}.eff_F);
                ESSlogLpCN(rep)   = summarypCN{i}.eff_LogL;
                TrainTimepCN(rep) = summarypCN{i}.elapsed;
                deltapCN(rep)     = summarypCN{i}.delta;
                accRatepCN(rep)   = summarypCN{i}.accRates.F;
            end
        end

        % --- pCNL ---
        fname = ['logGaussianCoxGirolami_pCNL_repeat' num2str(rep) '.mat'];
        if exist(fname,'file')
            load(fname);
            if iscell(summarypCNL) && length(summarypCNL) >= i
                ESSminpCNL(rep)    = min(summarypCNL{i}.eff_F);
                ESSmedianpCNL(rep) = median(summarypCNL{i}.eff_F);
                ESSmaxpCNL(rep)    = max(summarypCNL{i}.eff_F);
                ESSlogLpCNL(rep)   = summarypCNL{i}.eff_LogL;
                TrainTimepCNL(rep) = summarypCNL{i}.elapsed;
                deltapCNL(rep)     = summarypCNL{i}.delta;
                accRatepCNL(rep)   = summarypCNL{i}.accRates.F;
            end
        end

        % --- RHMC (N=64 only, i=1) ---
        if i == 1
            fname = ['logGaussianCoxGirolami_RHMC_repeat' num2str(rep) '_RHMC.mat'];
            if exist(fname,'file')
                load(fname);
                ESSminRHMC(rep)    = min(ess);
                ESSmedianRHMC(rep) = median(ess);
                ESSmaxRHMC(rep)    = max(ess);
                ESSlogLRHMC(rep)   = eff_LogL;
                TrainTimeRHMC(rep) = elapsedTime;
            end
        end

        % --- GI-MALA (N=64 only, i=1) ---
        if i == 1
            fname = ['logGaussianCoxGirolami_Marg_repeat' num2str(rep) '_gimala.mat'];
            if exist(fname,'file')
                load(fname);
                ESSminGimala(rep)    = min(ess);
                ESSmedianGimala(rep) = median(ess);
                ESSmaxGimala(rep)    = max(ess);
                ESSlogLGimala(rep)   = eff_LogL;
                TrainTimeGimala(rep) = elapsedTime;
                deltaGimala(rep)     = delta;
                accRateGimala(rep)   = accRates.F;
            end
        end

    end % rep loop

    % ---- LaTeX table ----
    fid = fopen([outdir tableNames{i} '.txt'],'w');
    fprintf(fid,'Method &  Time(s) & Acc. (\\%%)  & Step size $\\delta$  &  ESS (Min, Med, Max)  & Min ESS/s (1 st.d.) \\\\ \n');
    fprintf(fid,'\\midrule\n');

    if ~all(isnan(ESSminGimala))
        fprintf(fid,'GI-MALA  &  %1.1f  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimeGimala,'omitnan'), mean(accRateGimala,'omitnan'), mean(deltaGimala,'omitnan'), ...
            mean(ESSminGimala,'omitnan'), mean(ESSmedianGimala,'omitnan'), mean(ESSmaxGimala,'omitnan'), ...
            mean(ESSminGimala./TrainTimeGimala,'omitnan'), std(ESSminGimala./TrainTimeGimala,'omitnan'));
    end
    if ~all(isnan(ESSminMALA))
        fprintf(fid,'MALA  &  %1.1f  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimeMALA,'omitnan'), mean(accRateMALA,'omitnan'), mean(deltaMALA,'omitnan'), ...
            mean(ESSminMALA,'omitnan'), mean(ESSmedianMALA,'omitnan'), mean(ESSmaxMALA,'omitnan'), ...
            mean(ESSminMALA./TrainTimeMALA,'omitnan'), std(ESSminMALA./TrainTimeMALA,'omitnan'));
    end
    if ~all(isnan(ESSminMarg))
        fprintf(fid,'mGrad  &  %1.1f  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimeMarg,'omitnan'), mean(accRateMarg,'omitnan'), mean(deltaMarg,'omitnan'), ...
            mean(ESSminMarg,'omitnan'), mean(ESSmedianMarg,'omitnan'), mean(ESSmaxMarg,'omitnan'), ...
            mean(ESSminMarg./TrainTimeMarg,'omitnan'), std(ESSminMarg./TrainTimeMarg,'omitnan'));
    end
    if ~all(isnan(ESSminEllipt))
        fprintf(fid,'Ellipt  &  %1.1f  &   &   &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimeEllipt,'omitnan'), ...
            mean(ESSminEllipt,'omitnan'), mean(ESSmedianEllipt,'omitnan'), mean(ESSmaxEllipt,'omitnan'), ...
            mean(ESSminEllipt./TrainTimeEllipt,'omitnan'), std(ESSminEllipt./TrainTimeEllipt,'omitnan'));
    end
    if ~all(isnan(ESSminpCN))
        fprintf(fid,'pCN  &  %1.1f  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimepCN,'omitnan'), mean(accRatepCN,'omitnan'), mean(deltapCN,'omitnan'), ...
            mean(ESSminpCN,'omitnan'), mean(ESSmedianpCN,'omitnan'), mean(ESSmaxpCN,'omitnan'), ...
            mean(ESSminpCN./TrainTimepCN,'omitnan'), std(ESSminpCN./TrainTimepCN,'omitnan'));
    end
    if ~all(isnan(ESSminpCNL))
        fprintf(fid,'pCNL  &  %1.1f  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimepCNL,'omitnan'), mean(accRatepCNL,'omitnan'), mean(deltapCNL,'omitnan'), ...
            mean(ESSminpCNL,'omitnan'), mean(ESSmedianpCNL,'omitnan'), mean(ESSmaxpCNL,'omitnan'), ...
            mean(ESSminpCNL./TrainTimepCNL,'omitnan'), std(ESSminpCNL./TrainTimepCNL,'omitnan'));
    end
    if ~all(isnan(ESSminRHMC))
        fprintf(fid,'RHMC  &  %1.1f  &   &   &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimeRHMC,'omitnan'), ...
            mean(ESSminRHMC,'omitnan'), mean(ESSmedianRHMC,'omitnan'), mean(ESSmaxRHMC,'omitnan'), ...
            mean(ESSminRHMC./TrainTimeRHMC,'omitnan'), std(ESSminRHMC./TrainTimeRHMC,'omitnan'));
    end
    fclose(fid);
    fprintf('Written: %s%s.txt\n', outdir, tableNames{i});

end % i loop
