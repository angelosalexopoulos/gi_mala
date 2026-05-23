% make_results_binaryclassification.m
%
% Collects MCMC results for GP binary classification (fixed hyperparameters)
% and writes LaTeX table rows and log-likelihood trace figures.
%
% Methods included (Tables 3, 4, 9, 10, 11):
%   mGrad, pMALA, Ellipt, pCN, pCNL, GI-MALA
%
% Output written to diagrams/:
%   <Dataset>_table_fixedhypers.txt   LaTeX rows for Tables 3/4/9/10/11
%   <Dataset>_*_logL_fixedhypers.eps  Log-L traces for Figures 4/5
%
% Run from the repository root after run_all_experiments.m.

outdir = 'diagrams/';
fontsz = 26;
addpath results/;
addpath results/LogisticRegression_GP/;

Repeats = 10;

for datasetName = {'Australian' 'German' 'Heart' 'Pima' 'Ripley'}

    ds = datasetName{1};

    switch ds
        case 'Ripley',    ax = [1 10000 -90  -60 ];  tableNum = 11;  figNum = 0;
        case 'Pima',      ax = [1 10000 -260 -210];  tableNum = 10;  figNum = 0;
        case 'Heart',     ax = [1 10000 -120 -60 ];  tableNum = 3;   figNum = 5;
        case 'German',    ax = [1 10000 -500 -400];  tableNum = 9;   figNum = 0;
        case 'Australian',ax = [1 10000 -250 -150];  tableNum = 4;   figNum = 4;
    end

    % Pre-allocate all ESS arrays with NaN for robustness
    [ESSminMarg,  ESSmedianMarg,  ESSmaxMarg,  ESSlogLMarg,  TrainTimeMarg,  deltaMarg ] = deal(NaN(1,Repeats));
    [ESSminMALA,  ESSmedianMALA,  ESSmaxMALA,  ESSlogLMALA,  TrainTimeMALA,  deltaMALA ] = deal(NaN(1,Repeats));
    [ESSminEllipt,ESSmedianEllipt,ESSmaxEllipt,ESSlogLEllipt,TrainTimeEllipt           ] = deal(NaN(1,Repeats));
    [ESSminpCN,   ESSmedianpCN,   ESSmaxpCN,   ESSlogLpCN,   TrainTimepCN,   deltapCN  ] = deal(NaN(1,Repeats));
    [ESSminpCNL,  ESSmedianpCNL,  ESSmaxpCNL,  ESSlogLpCNL,  TrainTimepCNL,  deltapCNL ] = deal(NaN(1,Repeats));
    [ESSminGimala,ESSmedianGimala,ESSmaxGimala,ESSlogLGimala,TrainTimeGimala,deltaGimala] = deal(NaN(1,Repeats));

    for rep = 1:Repeats

        % --- mGrad ---
        fname = [ds '_repeat' num2str(rep) '_Marg_fixedhypers.mat'];
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
        fname = [ds '_repeat' num2str(rep) '_Marg_fixedhypers_pMALA.mat'];
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
        fname = [ds '_repeat' num2str(rep) '_Ellipt_fixedhypers.mat'];
        if exist(fname,'file')
            load(fname);
            ESSminEllipt(rep)    = min(summaryEllipt.eff_F);
            ESSmedianEllipt(rep) = median(summaryEllipt.eff_F);
            ESSmaxEllipt(rep)    = max(summaryEllipt.eff_F);
            ESSlogLEllipt(rep)   = summaryEllipt.eff_LogL;
            TrainTimeEllipt(rep) = summaryEllipt.elapsed;
        end

        % --- pCN ---
        fname = [ds '_repeat' num2str(rep) '_pCN_fixedhypers.mat'];
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
        fname = [ds '_repeat' num2str(rep) '_pCNL_fixedhypers.mat'];
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
        fname = [ds '_repeat' num2str(rep) '_Marg_fixedhypers_gimala.mat'];
        if exist(fname,'file')
            load(fname);
            ESSminGimala(rep)    = min(ess);
            ESSmedianGimala(rep) = median(ess);
            ESSmaxGimala(rep)    = max(ess);
            ESSlogLGimala(rep)   = eff_LogL;
            TrainTimeGimala(rep) = elapsedTime;
            deltaGimala(rep)     = delta;
        end

        % --- Log-likelihood trace figures (rep 1 only, Figure 4 = Australian, Figure 5 = Heart) ---
        if rep == 1 && figNum > 0
            pfx = [outdir 'Figure' num2str(figNum) '_'];

            if exist([ds '_repeat1_Marg_fixedhypers.mat'],'file')
                load([ds '_repeat1_Marg_fixedhypers.mat']);
                figure; plot(summaryMarg.LogL,'k');
                ylabel('Log-likelihood','fontsize',fontsz); axis(ax);
                set(gca,'fontsize',fontsz); title('mGrad');
                set(gca,'XTick',[3000 6000 9000]);
                exportgraphics(gcf,[pfx 'marg.pdf'],'Resolution',300);
            end

            if exist([ds '_repeat1_Marg_fixedhypers_pMALA.mat'],'file')
                load([ds '_repeat1_Marg_fixedhypers_pMALA.mat']);
                figure; plot(LogL,'k');
                xlabel('Sampling iteration','fontsize',fontsz);
                ylabel('Log-likelihood','fontsize',fontsz); axis(ax);
                set(gca,'fontsize',fontsz); title('pMALA');
                set(gca,'XTick',[3000 6000 9000]);
                exportgraphics(gcf,[pfx 'mala.pdf'],'Resolution',300);
            end

            if exist([ds '_repeat1_Ellipt_fixedhypers.mat'],'file')
                load([ds '_repeat1_Ellipt_fixedhypers.mat']);
                figure; plot(summaryEllipt.LogL,'k');
                xlabel('Sampling iteration','fontsize',fontsz); axis(ax);
                set(gca,'fontsize',fontsz); title('Ellipt');
                set(gca,'XTick',[3000 6000 9000]);
                exportgraphics(gcf,[pfx 'ellipt.pdf'],'Resolution',300);
            end

            if exist([ds '_repeat1_pCN_fixedhypers.mat'],'file')
                load([ds '_repeat1_pCN_fixedhypers.mat']);
                figure; plot(summarypCN.LogL,'k');
                xlabel('Sampling iteration','fontsize',fontsz);
                ylabel('Log-likelihood','fontsize',fontsz); axis(ax);
                set(gca,'fontsize',fontsz); title('pCN');
                set(gca,'XTick',[3000 6000 9000]);
                exportgraphics(gcf,[pfx 'pCN.pdf'],'Resolution',300);
            end

            if exist([ds '_repeat1_pCNL_fixedhypers.mat'],'file')
                load([ds '_repeat1_pCNL_fixedhypers.mat']);
                figure; plot(summarypCNL.LogL,'k');
                xlabel('Sampling iteration','fontsize',fontsz); axis(ax);
                set(gca,'fontsize',fontsz); title('pCNL');
                set(gca,'XTick',[3000 6000 9000]);
                exportgraphics(gcf,[pfx 'pCNL.pdf'],'Resolution',300);
            end

            if exist([ds '_repeat1_Marg_fixedhypers_gimala.mat'],'file')
                load([ds '_repeat1_Marg_fixedhypers_gimala.mat']);
                figure; plot(LogL,'k');
                xlabel('Sampling iteration','fontsize',fontsz);
                ylabel('Log-likelihood','fontsize',fontsz); axis(ax);
                set(gca,'fontsize',fontsz); title('GI-MALA');
                set(gca,'XTick',[3000 6000 9000]);
                exportgraphics(gcf,[pfx 'gimala.pdf'],'Resolution',300);
            end

        end % rep == 1 && figNum > 0
    end % rep loop

    % ---- LaTeX table ----
    fid = fopen([outdir 'Table' num2str(tableNum) '.txt'],'w');
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
    fprintf('Written: %sTable%d.txt  (%s)\n', outdir, tableNum, ds);

end % dataset loop
