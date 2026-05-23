% make_results_regression.m
%
% Collects MCMC results for GP regression with informative likelihood
% (regressinformlik_d1000.mat, sigma2 = 0.1^2, d ~ 1001) and writes
% the LaTeX table row for:
%
%   Table 15  -- mGrad, pMALA, Ellipt, pCN, pCNL, GI-MALA
%
% Columns: Time(s), step size delta, ESS (Min, Med, Max), Min ESS/s (1 st.d.)
%
% Output written to diagrams/:
%   Table15.txt
%
% Run from the repository root after run_all_experiments.m.

outdir = 'diagrams/';
addpath results/GPregression/;

Repeats = 10;

[ESSminMarg,   ESSmedianMarg,   ESSmaxMarg,   ESSlogLMarg,   TrainTimeMarg,   deltaMarg  ] = deal(NaN(1,Repeats));
[ESSminMALA,   ESSmedianMALA,   ESSmaxMALA,   ESSlogLMALA,   TrainTimeMALA,   deltaMALA  ] = deal(NaN(1,Repeats));
[ESSminEllipt, ESSmedianEllipt, ESSmaxEllipt, ESSlogLEllipt, TrainTimeEllipt             ] = deal(NaN(1,Repeats));
[ESSminpCN,    ESSmedianpCN,    ESSmaxpCN,    ESSlogLpCN,    TrainTimepCN,    deltapCN   ] = deal(NaN(1,Repeats));
[ESSminpCNL,   ESSmedianpCNL,   ESSmaxpCNL,   ESSlogLpCNL,   TrainTimepCNL,   deltapCNL  ] = deal(NaN(1,Repeats));
[ESSminGimala, ESSmedianGimala, ESSmaxGimala, ESSlogLGimala, TrainTimeGimala, deltaGimala] = deal(NaN(1,Repeats));

for rep = 1:Repeats

    % --- mGrad ---
    fname = ['regression_repeat' num2str(rep) '_Marg.mat'];
    if exist(fname,'file')
        load(fname);
        ESSminMarg(rep)    = min(summaryMarg.eff_F);
        ESSmedianMarg(rep) = median(summaryMarg.eff_F);
        ESSmaxMarg(rep)    = max(summaryMarg.eff_F);
        ESSlogLMarg(rep)   = summaryMarg.eff_LogL;
        TrainTimeMarg(rep) = summaryMarg.elapsed;
        deltaMarg(rep)     = summaryMarg.delta;
    end

    % --- pMALA ---
    fname = ['regression_repeat' num2str(rep) '_pMALA.mat'];
    if exist(fname,'file')
        load(fname);
        ESSminMALA(rep)    = min(ess);
        ESSmedianMALA(rep) = median(ess);
        ESSmaxMALA(rep)    = max(ess);
        ESSlogLMALA(rep)   = eff_LogL;
        TrainTimeMALA(rep) = elapsedTime;
        deltaMALA(rep)     = delta;
    end

    % --- Ellipt ---
    fname = ['regression_repeat' num2str(rep) '_Ellipt.mat'];
    if exist(fname,'file')
        load(fname);
        ESSminEllipt(rep)    = min(summaryEllipt.eff_F);
        ESSmedianEllipt(rep) = median(summaryEllipt.eff_F);
        ESSmaxEllipt(rep)    = max(summaryEllipt.eff_F);
        ESSlogLEllipt(rep)   = summaryEllipt.eff_LogL;
        TrainTimeEllipt(rep) = summaryEllipt.elapsed;
    end

    % --- pCN ---
    fname = ['regression_repeat' num2str(rep) '_pCN.mat'];
    if exist(fname,'file')
        load(fname);
        ESSminpCN(rep)    = min(summarypCN.eff_F);
        ESSmedianpCN(rep) = median(summarypCN.eff_F);
        ESSmaxpCN(rep)    = max(summarypCN.eff_F);
        ESSlogLpCN(rep)   = summarypCN.eff_LogL;
        TrainTimepCN(rep) = summarypCN.elapsed;
        deltapCN(rep)     = summarypCN.delta;
    end

    % --- pCNL ---
    fname = ['regression_repeat' num2str(rep) '_pCNL.mat'];
    if exist(fname,'file')
        load(fname);
        ESSminpCNL(rep)    = min(summarypCNL.eff_F);
        ESSmedianpCNL(rep) = median(summarypCNL.eff_F);
        ESSmaxpCNL(rep)    = max(summarypCNL.eff_F);
        ESSlogLpCNL(rep)   = summarypCNL.eff_LogL;
        TrainTimepCNL(rep) = summarypCNL.elapsed;
        deltapCNL(rep)     = summarypCNL.delta;
    end

    % --- GI-MALA ---
    fname = ['regression_repeat' num2str(rep) '_gimala.mat'];
    if exist(fname,'file')
        load(fname);
        ESSminGimala(rep)    = min(ess);
        ESSmedianGimala(rep) = median(ess);
        ESSmaxGimala(rep)    = max(ess);
        ESSlogLGimala(rep)   = eff_LogL;
        TrainTimeGimala(rep) = elapsedTime;
        deltaGimala(rep)     = delta;
    end

end

% ---- LaTeX table ----
fid = fopen([outdir 'Table15.txt'],'w');
fprintf(fid,' Method &  Time(s) & Step size $\\delta$  &  ESS (Min, Med, Max)  & Min ESS/s (1 st.d.) \\\\ \n');
fprintf(fid,'\\midrule\n');

if ~all(isnan(ESSminMarg))
    fprintf(fid,'mGrad  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
        mean(TrainTimeMarg,'omitnan'), mean(deltaMarg,'omitnan'), ...
        mean(ESSminMarg,'omitnan'), mean(ESSmedianMarg,'omitnan'), mean(ESSmaxMarg,'omitnan'), ...
        mean(ESSminMarg./TrainTimeMarg,'omitnan'), std(ESSminMarg./TrainTimeMarg,'omitnan'));
end
if ~all(isnan(ESSminMALA))
    fprintf(fid,'pMALA  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
        mean(TrainTimeMALA,'omitnan'), mean(deltaMALA,'omitnan'), ...
        mean(ESSminMALA,'omitnan'), mean(ESSmedianMALA,'omitnan'), mean(ESSmaxMALA,'omitnan'), ...
        mean(ESSminMALA./TrainTimeMALA,'omitnan'), std(ESSminMALA./TrainTimeMALA,'omitnan'));
end
if ~all(isnan(ESSminEllipt))
    fprintf(fid,'Ellipt  &  %1.1f  &   &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
        mean(TrainTimeEllipt,'omitnan'), ...
        mean(ESSminEllipt,'omitnan'), mean(ESSmedianEllipt,'omitnan'), mean(ESSmaxEllipt,'omitnan'), ...
        mean(ESSminEllipt./TrainTimeEllipt,'omitnan'), std(ESSminEllipt./TrainTimeEllipt,'omitnan'));
end
if ~all(isnan(ESSminpCN))
    fprintf(fid,'pCN  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
        mean(TrainTimepCN,'omitnan'), mean(deltapCN,'omitnan'), ...
        mean(ESSminpCN,'omitnan'), mean(ESSmedianpCN,'omitnan'), mean(ESSmaxpCN,'omitnan'), ...
        mean(ESSminpCN./TrainTimepCN,'omitnan'), std(ESSminpCN./TrainTimepCN,'omitnan'));
end
if ~all(isnan(ESSminpCNL))
    fprintf(fid,'pCNL  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
        mean(TrainTimepCNL,'omitnan'), mean(deltapCNL,'omitnan'), ...
        mean(ESSminpCNL,'omitnan'), mean(ESSmedianpCNL,'omitnan'), mean(ESSmaxpCNL,'omitnan'), ...
        mean(ESSminpCNL./TrainTimepCNL,'omitnan'), std(ESSminpCNL./TrainTimepCNL,'omitnan'));
end
fprintf(fid,'\\midrule\n');
if ~all(isnan(ESSminGimala))
    fprintf(fid,'GI-MALA  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
        mean(TrainTimeGimala,'omitnan'), mean(deltaGimala,'omitnan'), ...
        mean(ESSminGimala,'omitnan'), mean(ESSmedianGimala,'omitnan'), mean(ESSmaxGimala,'omitnan'), ...
        mean(ESSminGimala./TrainTimeGimala,'omitnan'), std(ESSminGimala./TrainTimeGimala,'omitnan'));
end
fclose(fid);
fprintf('Written: %sTable15.txt\n', outdir);
