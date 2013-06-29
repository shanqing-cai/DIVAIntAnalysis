function pert_analysis(subjID, varargin)
%% Constants
DATA_DIR = 'E:\DATA_CadLab\IntensityShift_normal';
L3DataDir = 'E:\DATA_CadLab\IntensityShift_normal\L3DATA';

gray = [0.5, 0.5, 0.5];
ALL_PHASES = {'Base', 'Ramp', 'Pert', 'Post'};

%% Options
bVis = isempty(fsic(varargin, '--noPlot'));

%%
check_dir(L3DataDir);
fns = dir(fullfile(L3DataDir, sprintf('%s (*).csv', subjID)));
assert(length(fns) == 1);
dfn = fullfile(L3DataDir, fns(1).name);
if ~isempty(strfind(dfn, '(DN)'))
    subjGrp = 'DN';
elseif ~isempty(strfind(dfn, '(UP)'))
    subjGrp = 'UP';
else
    error('Cannot determine group of subject');
end

fprintf(1, 'INFO: %s: Loading data from file: %s\n', mfilename, dfn);

dat = csvread(dfn);
assert(size(dat, 2) == 10);

%% 

% In the variable names here, "diff" means the difference between stressed
% and unstressed words
trialNum = dat(:, 1);
wordN = dat(:, 3);

WS_pI = dat(:, 5);
WU_pI = dat(:, 6);

phase = dat(:, 4);
fbDiff = dat(:, 9); 
micDiff = dat(:, 10);

diffDiff = fbDiff - micDiff;
diffPhase = diff(phase);
assert(isempty(find(diffPhase < 0)));

idxRamp0 = find(phase == 2, 1);
idxPert0 = find(phase == 3, 1);
idxPost0 = find(phase == 4, 1);

%% Load the sData and check proper trial correspondence between the csv and stat files
dDir = fullfile(DATA_DIR, subjGrp);
sData1 = intShift_readSubjData(fullfile(dDir, [subjID, '.stat']), ...
                              0, 'nStresses', 1, 'no_reverse');
sData1Rev = intShift_readSubjData(fullfile(dDir, [subjID, '.stat']), ...
                              0, 'nStresses', 1, 'reverse');

W1S_mI = nan(size(fbDiff));
W1U_mI = nan(size(fbDiff));

for i1 = 1 : numel(trialNum)
    if wordN(i1) ~= 1
        continue;
    end
    idx = find(sData1.trialNum == trialNum(i1));
    if isempty(idx)
        fprintf(2, 'WARNING: Cannot find trial number %d in sData from file %s\n', ...
                trialNum(i1), fullfile(dDir, [subjID, '.stat']));
        continue
    end
    
    assert(isequal(sData1.phase{idx}, ALL_PHASES{phase(i1)}));
    assert(isequal(sData1Rev.phase{idx}, ALL_PHASES{phase(i1)}));
    
    W1S_mI(i1) = sData1.mI(idx);
    W1U_mI(i1) = sData1Rev.mI(idx);
end

%% -- Reality check -- %
% figplot(W1S_mI - W1U_mI, micDiff, 'o');
% figplot(WS_pI - WU_pI, micDiff, 'o');
figplot(WS_pI(wordN == 1) - WU_pI(wordN == 1), micDiff(wordN == 1), 'mo');
hold on;
plot(WS_pI(wordN == 2) - WU_pI(wordN == 2), micDiff(wordN == 2), 'ro');
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
lims = [min([xs(1), ys(1)]), max([xs(2), ys(2)])];
plot(lims, lims, '-', 'Color', [0.5, 0.5, 0.5]);
set(gca, 'XLim', lims, 'YLim', lims);
axis square;
grid on;
xlabel('Stressed IL - unstressed IL (Column 5 - Column 6) (dB)');
ylabel('Column 10 (dB)');
legend({'W1 Stressed', 'W2 Stressed'});
title(strrep(subjID, '_', '\_'));
% figplot(WS_pI, W1S_mI);
% -- ~Reality check -- %

%% Visualization
if bVis
    % --- Values vs. trial number --- %
    figure('Position', [100, 100, 600, 600]);
    
    for i1 = 1 : 4
        subplot(4, 1, i1); hold on;
        
        if i1 == 1
            plot(trialNum, WS_pI, '-', 'Color', [0, 0, 0]);
            plot(trialNum, WU_pI, '-', 'Color', [0.5, 0.5, 0.5]);
            title('Peak intensity (mic)');
            legend({'S', 'U'}, 'FontSize', 6);
        elseif i1 == 2
            plot(trialNum, micDiff);
            title('Microphone: stressed - unstressed');
        elseif i1 == 3
            plot(trialNum, fbDiff);
            title('Feedback: stressed - unstressed');
        else
            plot(trialNum, diffDiff);
            title('Feedback (stressed - unstressed) - Mic (stressed - unstressed)');
        end
        set(gca, 'XLim', [0, trialNum(end) + 1]);

        xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
        if i1 == 1
            rng = [min([WS_pI; WU_pI]), max([WS_pI; WU_pI])];
            set(gca, 'YLim', [rng(1) - 0.05 * range(rng), rng(2) + 0.05 * range(rng)]);
            pause(0);
        else
            plot(xs, [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
        end

        ylabel('Diff. (dB)');

        plot(repmat(idxRamp0, 1, 2), ys, '-', 'Color', 'g');
        plot(repmat(idxPert0, 1, 2), ys, '-', 'Color', 'g');
        plot(repmat(idxPost0, 1, 2), ys, '-', 'Color', 'g');
    end
    
%     subplot(3, 1, 2); hold on;
%     plot(trialNum, fbDiff);
%     set(gca, 'XLim', [0, trialNum(end) + 1]);
%     
%     xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
%     plot(xs, [0, 0], '-', 'Color', gray);
%     
%     title('Feedback: stressed - unstressed');
%     xlabel('Trial number');
%     ylabel('Diff. (dB)');
%     
%     plot(repmat(idxRamp0, 1, 2), ys, '-', 'Color', 'g');
%     plot(repmat(idxPert0, 1, 2), ys, '-', 'Color', 'g');
%     plot(repmat(idxPost0, 1, 2), ys, '-', 'Color', 'g');
    
    % --- fbDiff vs. micDiff --- %
    figure();
    hold on;
    idxPert = find(phase == 3);
    plot(micDiff(idxPert), fbDiff(idxPert), 'o');
    xlabel('Microphone: stressed - unstressed (dB)');
    ylabel('Feedback: stressed - unstressed (dB)');
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    plot(xs, [0, 0], '-', 'Color', gray);
    plot([0, 0], ys, '-', 'Color', gray);
    box on;
    
    
    % --- fbDiff vs. diffDiff --- %
    figure('Position', [150, 150, 1200, 450]);    
    idxPert_W1 = find(phase == 3 & wordN == 1);
    idxPert_W2 = find(phase == 3 & wordN == 2);
    
    subplot(1, 3, 1);
    hold on;
    plot(micDiff(idxPert_W1), diffDiff(idxPert_W1), 'mo');
    plot(micDiff(idxPert_W2), diffDiff(idxPert_W2), 'ro');
    legend({'W1 stressed', 'W2 stressed'});
    
    xlabel('Microphone: stressed - unstressed (dB)');
    ylabel('Feedback (stressed - unstressed) - Mic (stressed - unstressed) (dB)');
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    plot(xs, [0, 0], '-', 'Color', gray);
    plot([0, 0], ys, '-', 'Color', gray);
    box on;
    
    subplot(1, 3, 2);
    hold on;
    plot(WU_pI(idxPert_W2), diffDiff(idxPert_W2), 'o');
    xlabel('Peak intensity of unstressed word in W2 trials (dB)');
    ylabel('Feedback (stressed - unstressed) - Mic (stressed - unstressed) (dB)');
    box on;
    
%     [sp_r, sp_t, sp_p] = spear(WU_pI(idxPert_W2), diffDiff(idxPert_W2));
    [lc_k, lc_r, lc_p] = lincorr(WU_pI(idxPert_W2), diffDiff(idxPert_W2));
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    plot(xs, lc_k(1) + lc_k(2) * xs, 'b-');
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.065 * range(ys), ...
         sprintf('Lin. corr. p = %.3e', lc_p));
    plot(xs, [0, 0], '-', 'Color', gray);
    set(gca, 'XLim', xs, 'YLim', ys);
    title('W2: Unstressed vs. pert. to contrast');
    
    subplot(1, 3, 3);
    hold on;
    plot(WS_pI(idxPert_W2), diffDiff(idxPert_W2), 'o');
    xlabel('Peak intensity of stressed word in W2 trials (dB)');
    ylabel('Feedback (stressed - unstressed) - Mic (stressed - unstressed) (dB)');
    box on;
    
%     [sp_r, sp_t, sp_p] = spear(WU_pI(idxPert_W2), diffDiff(idxPert_W2));
    [lc_k, lc_r, lc_p] = lincorr(WS_pI(idxPert_W2), diffDiff(idxPert_W2));
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    plot(xs, lc_k(1) + lc_k(2) * xs, 'b-');
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.065 * range(ys), ...
         sprintf('Lin. corr. p = %.3e', lc_p));
    plot(xs, [0, 0], '-', 'Color', gray);
    set(gca, 'XLim', xs, 'YLim', ys);
    title('W2: Stressed vs. pert. to contrast');
end


return
