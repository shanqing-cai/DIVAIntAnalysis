function sData = intShift_readSubjData(dataFN, bCD, varargin)
% txt = textread(dataFN, ...
%                '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s', ...
%                'delimiter', '\n');

DEBUG = 0;

nStresses = [1, 2];
if ~isempty(fsic(varargin, 'nStresses'))
    nStresses = varargin{fsic(varargin, 'nStresses') + 1};
end

bReverse = ~isempty(fsic(varargin, 'reverse'));

txt = textread(dataFN, ...
               '%s', ...
               'delimiter', '\n');
headers = regexp(txt{1}, '\t', 'split');

lenHeaders = length(headers);
idxEmptyHeaders = fsic(headers, '');
if DEBUG
    fprintf(1, 'INFO: length(headers) == %d\n', lenHeaders);
end

cTrialNum = strmatch('Trial', headers, 'exact');
cStress = strmatch('Stress', headers, 'exact');
cmF0_w1 = strmatch('W1_mF0', headers, 'exact');
cmI_w1 = strmatch('W1_mI', headers, 'exact');
cDur_w1 = strmatch('W1_dur', headers, 'exact');
cPhase = strmatch('Phase', headers, 'exact');
cStim = strmatch('Stimulus', headers, 'exact');
if isempty(cStim)
    cStim = strmatch('Stimuli', headers, 'exact');
end

cEpoch = strmatch('Epoch', headers, 'exact');

txt = txt(3 : end);

sData.trialNum = [];
sData.mF0 = [];     % Mean F0
sData.mI = [];  % Mean intensity
sData.dur = []; % Duration of stressed word
sData.nStress = [];
sData.phase = {};
sData.stim = [];
sData.epoch = [];

for i0 = 1 : numel(txt)
    tline = regexp(txt{i0}, '\t', 'split');
    if length(tline) ~= lenHeaders
        error('Length mismatch between tline and headers');
    end
    
    if ~isequal(fsic(tline, ''), idxEmptyHeaders)
        fprintf(2, 'WARNING: Mismatch in empty string locations between tline and headers\n');
    end
    
    
    nStress = tline{cStress};
    nStress = str2num(nStress(end));    
    nStim = tline{cStim};
    nStim = str2num(strrep(nStim, 'Stim', ''));
    nEpoch = tline{cEpoch};
    nEpoch = str2num(strrep(nEpoch, 'Epoch', ''));
    
    if isempty(find(nStresses == nStress))
        continue;
    end
    
    if nStress == 1
        nUnstress = 2;
    elseif nStress == 2
        nUnstress = 1;
    else
        error('Unexpected stress position: %d', nStress);
    end
    
    if nStress > 1
        pause(0);
    end
    
    sData.phase{end + 1} = tline{cPhase};
    sData.stim(end + 1) = nStim;
    sData.epoch(end + 1) = nEpoch;       
    sData.nStress(end + 1) = nStress;
    
    sData.trialNum(end + 1) = str2double(strrep(lower(tline(cTrialNum)), 'trial', ''));
    
    if bCD == 1
        sData.mF0(end + 1) = str2double(tline{cmF0_w1 + nStress - 1}) - ...
                             str2double(tline{cmF0_w1 + nUnstress - 1});
        sData.mI(end + 1) = str2double(tline{cmI_w1 + nStress - 1}) - ...
                            str2double(tline{cmI_w1 + nUnstress - 1});
        sData.dur(end + 1) = str2double(tline{cDur_w1 + nStress - 1}) - ...
                             str2double(tline{cDur_w1 + nUnstress - 1});
    else
        if ~bReverse
            sData.mF0(end + 1) = str2double(tline{cmF0_w1 + nStress - 1});
            sData.mI(end + 1) = str2double(tline{cmI_w1 + nStress - 1});
            sData.dur(end + 1) = str2double(tline{cDur_w1 + nStress - 1});
        else
            sData.mF0(end + 1) = str2double(tline{cmF0_w1 + nUnstress - 1});
            sData.mI(end + 1) = str2double(tline{cmI_w1 + nUnstress - 1});
            sData.dur(end + 1) = str2double(tline{cDur_w1 + nUnstress - 1});
        end
    end
end

if DEBUG
    fprintf(1, 'INFO: number of trials = %d\n', length(sData.phase));
end
return