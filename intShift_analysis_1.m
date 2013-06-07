function intShift_analysis_1(bCD, nStresses, varargin)
%% Optional input arguments
% cpa: 
%           Calculate composite prosody score.
% showByEpoch:
%           Show data by epochs
% showIndS: 
%           Show data from individual subjects
% reverse:
%           Reverse data between stressed (perturbed) and unstressed
%           (unpertrubed) words, e.g., for looking at the data from the
%           unstresed words
% permute N: 
%           Perform permutation test: N iterations (N > 0)

%% CONFIG
DATA_DIR = 'E:\DATA_CadLab\IntensityShift_normal';
% DATA_DIR = 'E:\DATA_CadLab\IntensityShift_normal\July2012_bak';

SDIRS = {'DN', 'UP'};
colors.DN = 'b';
colors.UP = 'r';
PHASES = {'Base', 'Ramp', 'Pert', 'Post'};
NEPOCHS = 60;

indS_figs_dir = 'E:\speechres\cadlab\indS_figs';

wPros = [1.0, 1.0, 1.0];    % For CPA: prosody weights: [w_intensity, w_F0 and w_dur]

P_UNC_THRESH_BY_EPOCH = 0.05;
P_PERMCORR_THRESH = 0.05;

%% Input arguments / options
bShowIndS = ~isempty(fsic(varargin, 'showIndS'));
bShowByEpoch = ~isempty(fsic(varargin, 'showByEpoch'));
bReverse = ~isempty(fsic(varargin, 'reverse'));
bCPA = ~isempty(fsic(varargin, 'cpa')) || ~isempty(fsic(varargin, 'CPA'));

if bCD && bReverse
    error('reverse can not be used with contrast distance');
end

nPerm = 0;
if ~isempty(fsic(varargin, 'permute'))
    nPerm = varargin{fsic(varargin, 'permute') + 1};
    assert(nPerm >= 1);
    
    check_dir('perm_files', '-create');
    
    permMatFN = fullfile('perm_files', sprintf('perm_dat_%d.mat', nPerm));
end

%% Load data

sIDs = struct;
mF0_byP = struct; % Mean F0, by phase
mI_byP = struct; % Mean intensity, by phase
dur_byP = struct; % duration, by phsae

rmF0_byP = struct; % Relative mean F0, by phase
rmI_byP = struct; % Relative mean intensity, by phase
rdur_byP = struct; % Relative duration, by phase

mF0_byE = struct; % Mean F0, by phase
mI_byE = struct; % Mean intensity, by phase
dur_byE = struct; % duration, by phsae

rmF0_byE = struct;
rmI_byE = struct;
rdur_byE = struct;

for i0 = 1 : numel(SDIRS)    
    sDir = SDIRS{i0};
    
    sIDs.(sDir) = {};
    mF0_byP.(sDir) = nan(0, 4); % Four phases: Base, Ramp, Pert and Post
    mI_byP.(sDir) = nan(0, 4);
    dur_byP.(sDir) = nan(0, 4);
    
    rmF0_byP.(sDir) = nan(0, 4); % Four phases: Base, Ramp, Pert and Post
    rmI_byP.(sDir) = nan(0, 4);
    rdur_byP.(sDir) = nan(0, 4);
    
    mF0_byE.(sDir) = nan(0, NEPOCHS); % Four phases: Base, Ramp, Pert and Post
    mI_byE.(sDir) = nan(0, NEPOCHS);
    dur_byE.(sDir) = nan(0, NEPOCHS);
    
    rmF0_byE.(sDir) = nan(0, NEPOCHS);
    rmI_byE.(sDir) = nan(0, NEPOCHS);
    rdur_byE.(sDir) = nan(0, NEPOCHS);
    
    dDir = fullfile(DATA_DIR, sDir);
    
    stat_fns = dir(fullfile(dDir, '*.stat'));
    
    for i1 = 1 : numel(stat_fns)
        sID = strrep(stat_fns(i1).name, '.stat', '');
        sIDs.(sDir){end + 1} = sID;
        fprintf('Loading data from subject %s: %s...\n', sDir, sID);
        
        if bReverse
            revWord = 'reverse';
        else
            revWord = 'no_reverse';
        end
        
        sData = intShift_readSubjData(fullfile(dDir, stat_fns(i1).name), ...
                                      bCD, 'nStresses', nStresses, revWord);
        
        idx_base = strmatch('Base', sData.phase, 'exact');
        idx_ramp = strmatch('Ramp', sData.phase, 'exact');
        idx_pert = strmatch('Pert', sData.phase, 'exact');
        idx_post = strmatch('Post', sData.phase, 'exact');
        
        % By-phase data
        p_mF0 = [nanmean(sData.mF0(idx_base)), ...
                 nanmean(sData.mF0(idx_ramp)), ...
                 nanmean(sData.mF0(idx_pert)), ...
                 nanmean(sData.mF0(idx_post))];
        p_mI = [nanmean(sData.mI(idx_base)), ...
                nanmean(sData.mI(idx_ramp)), ...
                nanmean(sData.mI(idx_pert)), ...
                nanmean(sData.mI(idx_post))];
        p_dur = [nanmean(sData.dur(idx_base)), ...
                nanmean(sData.dur(idx_ramp)), ...
                nanmean(sData.dur(idx_pert)), ...
                nanmean(sData.dur(idx_post))];
        
        mF0_byP.(sDir) = [mF0_byP.(sDir); p_mF0];
        mI_byP.(sDir) = [mI_byP.(sDir); p_mI];
        dur_byP.(sDir) = [dur_byP.(sDir); p_dur];
        
        % Compute relative changes from the Base phase (by phase)
        rmF0_byP.(sDir) = [rmF0_byP.(sDir); p_mF0 / p_mF0(1)];
        rmI_byP.(sDir) = [rmI_byP.(sDir); p_mI / p_mI(1)];
        rdur_byP.(sDir) = [rdur_byP.(sDir); p_dur / p_dur(1)];
        
        % By-epoch data
        e_mF0 = nan(1, NEPOCHS);
        e_mI = nan(1, NEPOCHS);
        e_dur = nan(1, NEPOCHS);
        for i2 = 1 : NEPOCHS
            idxe = find(sData.epoch == i2);
            e_mF0(i2) = nanmean(sData.mF0(idxe));
            e_mI(i2) = nanmean(sData.mI(idxe));
            e_dur(i2) = nanmean(sData.dur(idxe));
        end
        
        mF0_byE.(sDir) = [mF0_byE.(sDir); e_mF0];
        mI_byE.(sDir) = [mI_byE.(sDir); e_mI];
        dur_byE.(sDir) = [dur_byE.(sDir); e_dur];
        
        % Compute relative change from the base phase (by epoch)
        rmF0_byE.(sDir) = [rmF0_byE.(sDir); e_mF0 / p_mF0(1)];
        rmI_byE.(sDir) = [rmI_byE.(sDir); e_mI / p_mI(1)];
        rdur_byE.(sDir) = [rdur_byE.(sDir); e_dur / p_dur(1)];
    end
end

%% Optional:  Calculate the composite prosody adaptation (CPA) score
if bCPA
    cpa_byP = struct;
    cpa_byE = struct;
    for i1 = 1 : numel(SDIRS)
        sDir = SDIRS{i1};

        cpa_byP.(sDir) = (rmI_byP.(sDir) * wPros(1) + rmF0_byP.(sDir) * wPros(2) + rdur_byP.(sDir) * wPros(3)) / sum(wPros);
        cpa_byE.(sDir) = (rmI_byE.(sDir) * wPros(1) + rmF0_byE.(sDir) * wPros(2) + rdur_byE.(sDir) * wPros(3)) / sum(wPros);

    end
end

%% 
nPhases = numel(PHASES);
if bShowIndS
    figure('Position', [50, 50, 1200, 350]);
    for j0 = 1 : 3 
        if j0 == 1;         meas = rmI_byP;     measName = 'intensity';
        elseif j0 == 2;     meas = rmF0_byP;    measName = 'F0';
        elseif j0 == 3;     meas = rdur_byP;    measName = 'duration';
        end
        
        subplot(1, 3, j0);
        hold on;
        plot([0, nPhases + 1], [1, 1], '-', ...
             'Color', [0.5, 0.5, 0.5]);
        for i0 = 1 : numel(SDIRS)
            sDir = SDIRS{i0};
            plot(meas.(sDir)', 'o-', 'Color', colors.(sDir));
        end
        
        xlabel('Phase');
        if bCD
            ylabel(['Normalized ', measName, ' C.D.']);
        else
            ylabel(['Normalized ', measName]);
        end
        
        set(gca, 'XLim', [0, nPhases + 1]);
        set(gca, 'XTick', [1 : nPhases]);
        set(gca, 'XTickLabel', PHASES);
    end
end



%% Visualization and some comparisons
figure('Position', [50, 250, 1200, 350]);

nTot = numel(sIDs.DN) + numel(sIDs.UP);
if nPerm > 0 && isfile(permMatFN)
    load(permMatFN); 
    if exist('rp_ps', 'var') && exist('rp_ts', 'var') ...
       && size(rp_ps, 1) == nPerm && size(rp_ts, 1) == nPerm
        bPerm = 0;
    else
        bPerm = 1;
    end
else
    bPerm = 1;
end
if bPerm == 1
    rp_ps = nan(nPerm, 4, 3 + bCPA); % Permutations, phases, # of tests
    rp_ts = nan(nPerm, 4, 3 + bCPA);
end

unc_ps = nan(4, 3 + bCPA); % Phases, # of tests

for k0 = 0 : 1 : nPerm * bPerm
    % --- Random permutation --- %
    if k0 >= 1
        if k0 == 1
            nc = print_progress_bar(0, nPerm, sprintf('Performing permutation test (by phase)'));
        end
        rpidx = randperm(nTot);
    end
    
    for j0 = 1 : 3 + bCPA
        if k0 == 0
            if j0 == 4
                figure('Name', 'Composite prosody adaptation');            
            else           
                subplot(1, 3, j0);
            end
        end

        hold on;

        if j0 == 1
            meas = rmI_byP;
            meas1 = mI_byP;
            measName = 'intensity';
        elseif j0 == 2
            meas = rmF0_byP;
            meas1 = mF0_byP;
            measName = 'F0';
        elseif j0 == 3
            meas = rdur_byP;
            meas1 = dur_byP;
            measName = 'duration';
        elseif j0 == 4
            meas = cpa_byP;
            meas1 = [];
            measName = 'CPA';
        end
        
        % --- Random permutation --- %
        if k0 >= 1
            meas_a = [meas.DN; meas.UP];
            meas_a = meas_a(rpidx, :);
            meas.DN = meas_a(1 : size(meas.DN, 1), :);
            meas.UP = meas_a(size(meas.DN, 1) + 1 : nTot, :);
        end

        if k0 == 0
            for i0 = 1 : numel(SDIRS)
                sDir = SDIRS{i0};
                errorbar(1 : size(meas.(sDir), 2), ...
                         mean(meas.(sDir)), ste(meas.(sDir)), ...
                         'o-', 'Color', colors.(sDir));
            end
            xlabel('Phase');
            if bCD
                ylabel(['Normalized ', measName, ' C.D.']);
            else
                ylabel(['Normalized ', measName]);
            end
            title(measName);

            set(gca, 'XLim', [0, nPhases + 1]);
            set(gca, 'XTick', [1 : nPhases]);
            set(gca, 'XTickLabel', PHASES);
            
            ys = get(gca, 'YLim');
        end

        % Between-group, same-phase t-tests        
        for i1 = 2 : numel(PHASES)
            [~, p_t2, ~, t2_stats] = ttest2(meas.(SDIRS{1})(:, i1), meas.(SDIRS{2})(:, i1));
            
            if k0 == 0
                text(i1 - 0.3, ys(2) - 0.05 * range(ys), sprintf('p=%.3f', p_t2));
                
                unc_ps(i1, j0) = p_t2;
            else
                rp_ps(k0, i1, j0) = p_t2;
                rp_ts(k0, i1, j0) = t2_stats.tstat;
                
                if mod(k0, round(nPerm / 10)) == 0
                    for k1 = 1 : nc; fprintf(1, '\b'); end 
                    print_progress_bar(k0, nPerm, sprintf('Performing permutation test (by phase)'));
                end
            end
        end
        
        if k0 == 0
            text(0.1, ys(2) - 0.05 * range(ys), 'b/w group:');

            legend(SDIRS, 'Location', 'Southwest');
            plot([0, nPhases + 1], [1, 1], '-', 'Color', [0.5, 0.5, 0.5]);        

            % Within-group, between-phase t-tests
            for i0 = 1 : numel(SDIRS)
                sDir = SDIRS{i0};
                ps_t = nan(1, numel(PHASES));
                mean_meas = mean(meas.(sDir));  
                ys = get(gca, 'YLim');
                for i1 = 2 : numel(PHASES)
                    if j0 < 4
                        [h_foo, ps_t(i1)] = ttest(meas1.(sDir)(:, i1), meas1.(sDir)(:, 1));
                    else
                        [h_foo, ps_t(i1)] = ttest(meas.(sDir)(:, i1) - 1);
                    end

                    if ps_t(i1) < 0.05
                        fw = 'bold';
                    else
                        fw = 'light';
                    end
                    text(i1 + 0.07, mean_meas(i1) - 0.02 * range(ys), ...
                         sprintf('%.3f', ps_t(i1)), ...
                         'Color', colors.(sDir), 'FontSize', 9, 'FontWeight', fw);
                end        
            end

            % Within-group, RM-ANOVA with post-hoc Tukey
            if j0 < 4
                for i1 = 1 : numel(SDIRS)
                    sDir = SDIRS{i1};
                    rma_res = RM_ANOVA_1W(meas1.(sDir), PHASES, ...
                                          'contrasts', {[-1, 1, 0, 0], ...
                                                        [-1, 0, 1, 0], ...
                                                        [-1, 0, 0, 1]});
                    fprintf('%s: sDir = %s: RM-ANOVA results:\n', ...
                            measName, sDir);
                    fprintf('Omnibus: F(%d, %d) = %.3f, p = %3f\n', ...
                            rma_res.omniRes.df_A, rma_res.omniRes.df_SA, ...
                            rma_res.omniRes.F, rma_res.omniRes.p);
                    fprintf('Post-hoc Tukey HSD: \n');
                    fprintf('\tBase vs. Ramp: h = %d\n', ...
                            rma_res.tukeyRes{1}.h);
                    fprintf('\tBase vs. Pert: h = %d\n', ...
                            rma_res.tukeyRes{2}.h);
                    fprintf('\tBase vs. Post: h = %d\n', ...
                            rma_res.tukeyRes{3}.h);
                end
                fprintf('\n');
            end 
        end
    end
end

if nPerm > 0 && bPerm
    fprintf(1, '\n');
    
    save(permMatFN, 'rp_ps', 'rp_ts');
    check_file(permMatFN);
    fprintf(1, 'INFO: random permutation data saved to file: %s\n\n', permMatFN);
end

% --- Print permutation test results (by phase, between-group) --- %
if nPerm > 0
    itemNames = {'Intensity', 'F0', 'duration', 'CPA'};
    
    fprintf(1, '=== Permutation-corrected p-values, by-phase: ===\n');
    
    for i1 = 1 : 3 + bCPA        
        t_ps = rp_ps(:, 2 : end, i1);
        min_t_ps = min(t_ps');
        fprintf(1, '\t== %s: ==\n', itemNames{i1});
        
        for i2 = 2 : numel(PHASES)           
            t_corr_p = numel(find(t_ps <= unc_ps(i2, i1))) / nPerm;
            fprintf(1, '\t\t%s: p < %.4f', PHASES{i2}, t_corr_p);
            if t_corr_p < P_PERMCORR_THRESH
                fprintf(1, ' *');
            end
            fprintf(1, '\n');
        end
    end
end

%% Visualization: by Epoch
if bShowByEpoch
    figure('Position', [50, 50, 1400, 350]);

    for j0 = 1 : 3 + bCPA
        if j0 == 4
            figure('Name', 'CPA by epoch');
        else
            subplot(1, 3, j0);
        end
        hold on;

        if j0 == 1
            meas = rmI_byE;
            meas1 = mI_byE;
            measName = 'intensity';
        elseif j0 == 2
            meas = rmF0_byE;
            meas1 = mF0_byE;
            measName = 'F0';
        elseif j0 == 3
            meas = rdur_byE;
            meas1 = dur_byE;
            measName = 'duration';
        else
            meas = cpa_byE;
            meas1 = [];
            measName = 'CPA';
        end

        for i0 = 1 : numel(SDIRS)
            sDir = SDIRS{i0};
            errorbar(1 : size(meas.(sDir), 2), ...
                     mean(meas.(sDir)), ste(meas.(sDir)), ...
                     '-', 'Color', colors.(sDir));
        end
        xlabel('Phase');
        if bCD
            ylabel(['Normalized ', measName, ' C.D.']);
        else
            ylabel(['Normalized ', measName]);
        end
        title(measName);

        set(gca, 'XLim', [0, NEPOCHS + 1]);
%         set(gca, 'XTick', [1 : NEPOCHS]);        

        % Between-group, same-phase t-tests
        ys = get(gca, 'YLim');
        for i1 = 2 : NEPOCHS
            % -- Between-group significance -- %
            [~, p_t2] = ttest2(meas.(SDIRS{1})(:, i1), meas.(SDIRS{2})(:, i1));
            if p_t2 < P_UNC_THRESH_BY_EPOCH
                plot(i1, ys(2) - 0.04 * range(ys), 'k*', 'MarkerSize', 7);
            end
            
            % -- Within-group significance -- %
            for i2 = 1 : numel(SDIRS)
                sDir = SDIRS{i2};
                
                [~, p_t1] = ttest(meas.(sDir)(:, i1) - 1);
                if p_t1 < P_UNC_THRESH_BY_EPOCH
                    mkFaceClr = colors.(sDir);
                else
                    mkFaceClr = [1, 1, 1];
                end
                
                if j0 < 4
                    mkSize = 5;
                else
                    mkSize = 7;
                end
                
                plot(i1, mean(meas.(sDir)(:, i1)), 'o', ...
                     'MarkerEdgeColor', colors.(sDir), ...
                     'MarkerFaceColor', mkFaceClr, ...
                     'MarkerSize', mkSize);
            end
%             text(i1 - 0.3, ys(2) - 0.05 * range(ys), sprintf('%.2f', p_t2), 'FontSize', 6);
        end
        text(0.1, ys(2) - 0.04 * range(ys), 'b/w group:');

        legend(SDIRS, 'Location', 'Southwest');
        plot([0, NEPOCHS + 1], [1, 1], '-', 'Color', [0.5, 0.5, 0.5]);

        % Within-group, between-phase t-tests
%         for i0 = 1 : numel(SDIRS)
%             sDir = SDIRS{i0};
%             ps_t = nan(1, numel(PHASES));
%             mean_meas = mean(meas.(sDir));  
%             ys = get(gca, 'YLim');
%             for i1 = 2 : numel(PHASES)            
%                 [h_foo, ps_t(i1)] = ttest(meas1.(sDir)(:, i1), meas1.(sDir)(:, 1));
% 
%                 if ps_t(i1) < 0.05
%                     fw = 'bold';
%                 else
%                     fw = 'light';
%                 end
%                 text(i1 + 0.07, mean_meas(i1) - 0.02 * range(ys), ...
%                      sprintf('%.3f', ps_t(i1)), ...
%                      'Color', colors.(sDir), 'FontSize', 9, 'FontWeight', fw);
%             end        
%         end

        % Within-group, RM-ANOVA with post-hoc Tukey
%         for i1 = 1 : numel(SDIRS)
%             sDir = SDIRS{i1};
%             rma_res = RM_ANOVA_1W(meas1.(sDir), PHASES, ...
%                                   'contrasts', {[-1, 1, 0, 0], ...
%                                                 [-1, 0, 1, 0], ...
%                                                 [-1, 0, 0, 1]});
%             fprintf('%s: sDir = %s: RM-ANOVA results:\n', ...
%                     measName, sDir);
%             fprintf('Omnibus: F(%d, %d) = %.3f, p = %3f\n', ...
%                     rma_res.omniRes.df_A, rma_res.omniRes.df_SA, ...
%                     rma_res.omniRes.F, rma_res.omniRes.p);
%             fprintf('Post-hoc Tukey HSD: \n');
%             fprintf('\tBase vs. Ramp: h = %d\n', ...
%                     rma_res.tukeyRes{1}.h);
%             fprintf('\tBase vs. Pert: h = %d\n', ...
%                     rma_res.tukeyRes{2}.h);
%             fprintf('\tBase vs. Post: h = %d\n', ...
%                     rma_res.tukeyRes{3}.h);
%         end
%         fprintf('\n');

    end
    
    if bShowIndS        
        groups = fields(mI_byE);
        
        for k1 = 1 : numel(groups)
            grp = groups{k1};
                        
            for k2 = 1 : size(mI_byE.(grp), 1)
                t_sID = sIDs.(grp){k2};                
                hFigEpochIndS = figure('Position', [80, 80, 1200, 400], ...
                                       'Name', sprintf('IndS byE: %s', t_sID));
                                  
                for k3 = 1 : 3
                    if k3 == 1
                        iem = mI_byE.(grp)(k2, :);
                        measName = 'Epoch-mean intensity';
                    elseif k3 == 2
                        iem = mF0_byE.(grp)(k2, :);
                        measName = 'Epoch-mean F0 (Hz)';
                    else
                        iem = dur_byE.(grp)(k2, :);
                        measName = 'Epoch-mean duration (ms)';
                    end
                    
                    subplot(1, 3, k3); 
                    axis tight;
                    plot(1 : length(iem), iem, 'bo-'); 
                    ylabel(measName);
                    xlabel('Epoch #');
                    hold on;
                    
                    if k3 == 2
                        title(sprintf('%s (%s)', strrep(t_sID, '_', '\_'), grp));
                    end
                end                                
                                
                figFN = fullfile(indS_figs_dir, [t_sID, '.tif']);
                if isfile(figFN)
                    delete(figFN);
                end
                saveas(hFigEpochIndS, figFN, 'tif');
                if isfile(figFN)
                    fprintf(1, 'Saved the plot of subject %s to file: %s\n', t_sID, figFN);
                else
                    error('Failed to save to image file: %s (%s)', figFN, grp);
                end
                
                delete(hFigEpochIndS);
                drawnow;
            end                        
        end
    end
end

%% Write to excel file (for SYSTAT use);
for j0 = 1 : 3
    if j0 == 1
        meas = mI_byP;
        measName = 'intensity';
    elseif j0 == 2
        meas = mF0_byP;
        measName = 'F0';
    elseif j0 == 3
        meas = dur_byP;
        measName = 'duration';
    end
    
    dataMat = [[1 * ones(numel(sIDs.(SDIRS{1})), 1), meas.(SDIRS{1})]; ...
               [2 * ones(numel(sIDs.(SDIRS{2})), 1), meas.(SDIRS{2})]];
    dataCell = [{'SDIR'}, PHASES];
    dataCell = [dataCell; num2cell(dataMat)];
    
    if bCD
        xls_fn = fullfile(DATA_DIR, sprintf('%s_cd.xls', measName));
    else
        xls_fn = fullfile(DATA_DIR, sprintf('%s.xls', measName));
    end
    xlswrite(xls_fn, dataCell);
    fprintf('%s data written to file %s\n', measName, xls_fn);
end
return