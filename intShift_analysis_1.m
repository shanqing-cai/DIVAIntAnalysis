function intShift_analysis_1(bCD, nStresses, varargin)
%% CONFIG
DATA_DIR = 'E:\DATA_CadLab\IntensityShift_normal';
% DATA_DIR = 'E:\DATA_CadLab\IntensityShift_normal\July2012_bak';

SDIRS = {'DN', 'UP'};
colors.DN = 'b';
colors.UP = 'r';
PHASES = {'Base', 'Ramp', 'Pert', 'Post'};
NEPOCHS = 60;

indS_figs_dir = 'E:\speechres\cadlab\indS_figs';

%% Input arguments / options
bShowIndS = ~isempty(fsic(varargin, 'showIndS'));
bShowByEpoch = ~isempty(fsic(varargin, 'showByEpoch'));
bReverse = ~isempty(fsic(varargin, 'reverse'));

if bCD && bReverse
    error('reverse can not be used with contrast distance');
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
figure('Position', [50, 50, 1200, 350]);

for j0 = 1 : 3
    subplot(1, 3, j0);
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
    end
    
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
    
    % Between-group, same-phase t-tests
    ys = get(gca, 'YLim');
    for i1 = 2 : numel(PHASES)
        [h_t2, p_t2] = ttest2(meas.(SDIRS{1})(:, i1), meas.(SDIRS{2})(:, i1));        
        text(i1 - 0.3, ys(2) - 0.05 * range(ys), sprintf('p=%.3f', p_t2));
    end
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
            [h_foo, ps_t(i1)] = ttest(meas1.(sDir)(:, i1), meas1.(sDir)(:, 1));
            
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

%% Visualization: by Epoch
if bShowByEpoch
    figure('Position', [50, 50, 1400, 350]);

    for j0 = 1 : 3
        subplot(1, 3, j0);
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
        end

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

        set(gca, 'XLim', [0, NEPOCHS + 1]);
%         set(gca, 'XTick', [1 : NEPOCHS]);        

        % Between-group, same-phase t-tests
        ys = get(gca, 'YLim');
        for i1 = 2 : NEPOCHS
            [h_t2, p_t2] = ttest2(meas.(SDIRS{1})(:, i1), meas.(SDIRS{2})(:, i1));        
%             text(i1 - 0.3, ys(2) - 0.05 * range(ys), sprintf('%.2f', p_t2), 'FontSize', 6);
        end
        text(0.1, ys(2) - 0.05 * range(ys), 'b/w group:');

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