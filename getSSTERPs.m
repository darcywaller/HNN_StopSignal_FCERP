%% Analysis code associated with "Biophysical microcircuit mechanisms of 
% frontocentral evoked potentials support a race model interpretation of
% the Stop-Signal P3"

% Written by Darcy Diesburg (2023), Brown University, with P3 onset
% analysis replicated from script by Jan Wessel (i.e., Wessel and Aron, 2015).
% Contact: darcy.diesburg@gmail.com

%% Generate P3 ERPs from ALLSST datasets
% find all files
clear; clc; close all;
folders.data = '/gpfs/data/sjones/shared/StopSignal/EEG data/'; % set to location of OSF data from Wessel 2020
folders.base = pwd; 
folders.stats = fullfile(folders.base,'/stats'); mkdir(folders.stats);
addpath(folders.base,'functions');

% Requires EEGLab - set path to it
addpath('/gpfs/home/ddiesbur/eeglab/'); eeglab;

% get datasets
fname = dir(fullfile(folders.data,'*.set'));
fname = {fname.name};
% exclusion (for participants to be excluded after auto-preprocessing)
exclude = 70; % 70 = Subject161.ft. Didn't follow task instructions (Never stopped).
fname(exclude) = [];

%% make ERPs
% preallocate
win = [-.5 .5];
FAIL_STOPS = zeros(length(fname),500);
SUCC_STOPS = FAIL_STOPS;

for is = 1:length(fname)
    clear EEG bestcomp

    % load in EEG sets pp'ed from the ALLSST folder
    EEG = pop_loadset(fullfile(folders.data,fname{is}));
    substring = strrep(strjoin(regexp(fname{is},'[0-9]','match')),' ','');
    
    % run algorithmic detection of FC-P3 onset comp, correct sets where the
    % component isn't clear enough for algorithmic detection

    if strcmpi(substring,'10')
        bestcomp = 11;
    elseif strcmpi(substring,'100')
        bestcomp = 3;
    elseif strcmpi(substring,'128') || strcmpi(substring,'233')
        bestcomp = 6;
    elseif strcmpi(substring,'130')
        bestcomp = 13;
    elseif strcmpi(substring,'132') || strcmpi(substring,'16') || strcmpi(substring,'166') || strcmpi(substring,'91') || strcmpi(substring,'77')
        bestcomp = 1;
    elseif strcmpi(substring,'153') || strcmpi(substring,'188') || strcmpi(substring,'224') || strcmpi(substring,'48')
        bestcomp = 5;
    elseif strcmpi(substring,'169') || strcmpi(substring,'192') || strcmpi(substring,'90')
        bestcomp = 2;
    elseif strcmpi(substring,'186')
        bestcomp = 23;
    elseif strcmpi(substring,'189') || strcmpi(substring,'199') || strcmpi(substring,'202') || strcmpi(substring,'22')
        bestcomp = 12;
    elseif strcmpi(substring,'223')
        bestcomp = 10;
    elseif strcmpi(substring,'37')
        bestcomp = 14;
    elseif strcmpi(substring,'54')
        bestcomp = 16;
    elseif strcmpi(substring,'56')
        bestcomp = 4;
    elseif strcmpi(substring,'57')
        bestcomp = 9;
    elseif strcmpi(substring,'96') || strcmpi(substring,'144')
        bestcomp = 7;
    elseif strcmpi(substring,'92')
        bestcomp = 8;
    else
        bestcomp = cipi(EEG);
    end
    % get P3 IC comp
    icweights = 1:size(EEG.icaweights,1);
    
    % extract time series for that comp specifically
    EEG = pop_select(EEG,'nochannel',{'HEOG','VEOG'});    
    EEG = pop_subcomp(EEG, setxor(icweights,bestcomp));

    % epoch
    EEG = pop_epoch(EEG,{'S200'},win);
    EEG = pop_rmbase(EEG, [-100 0]);
    sEEG = pop_selectevent(EEG,'type',{'S200'},'acc',3,'deleteevents','on','deleteepochs','on','invertepochs','off');
    fEEG = pop_selectevent(EEG,'type',{'S200'},'acc',4,'deleteevents','on','deleteepochs','on','invertepochs','off');
    
    % save channel data
    succ_avg_FCZ = mean(mean(sEEG.data(find(strcmpi({EEG.chanlocs.labels},'FCz')),:,:),1),3);
    fail_avg_FCZ = mean(mean(fEEG.data(find(strcmpi({EEG.chanlocs.labels},'FCz')),:,:),1),3);
    
    SUCC_STOPS(is,:) = succ_avg_FCZ; 
    FAIL_STOPS(is,:) = fail_avg_FCZ;

end

%% Make ERP data into text files for HNN modeling software
% each line is a sampling point with voltage value at given timepoint
time_vector = 1:2:1000;
volt_vector = mean(SUCC_STOPS(:,251:end),1); % just mean erp starting at 0
fail_volt_vector = mean(FAIL_STOPS(:,251:end),1);
neg_volt_vector = -mean(SUCC_STOPS(:,251:end),1);

succ_fileID = fopen('SUCCSTOP.txt','w');
for il = 1:length(volt_vector)
    if il == 1
        string2print = [num2str(time_vector(il)) ' ' num2str(volt_vector(il))];
    else
        string2print = ['\n' num2str(time_vector(il)) ' ' num2str(volt_vector(il))];
    end
    fprintf(succ_fileID,string2print);
end
fclose(succ_fileID);

fail_fileID = fopen('FAILSTOP.txt','w');
for il = 1:length(fail_volt_vector)
    if il == 1
        string2print = [num2str(time_vector(il)) ' ' num2str(fail_volt_vector(il))];
    else
        string2print = ['\n' num2str(time_vector(il)) ' ' num2str(fail_volt_vector(il))];
    end
    fprintf(fail_fileID,string2print);
end
fclose(fail_fileID);

succ_fileID = fopen('SUCCSTOP_neg.txt','w');
for il = 1:length(neg_volt_vector)
    if il == 1
        string2print = [num2str(time_vector(il)) ' ' num2str(neg_volt_vector(il))];
    else
        string2print = ['\n' num2str(time_vector(il)) ' ' num2str(neg_volt_vector(il))];
    end
    fprintf(succ_fileID,string2print);
end
fclose(succ_fileID);
%% In the following sections, statistically identify P3 onset
% We use the procedure described in Wessel and Aron 2015

%% Grab ERPs for SS and FS as well as matched go trials with similar SSD
baseline = [-100 0];
plotlimits = [-100 1000];
mc = 1000;
trialnumbers = zeros(length(fname),3);
mcmc = 1; % bootstrapping on or off

for is = 1:length(fname)
    
    % load EEG file
    EEG = pop_loadset(fullfile(folders.data,fname{is}));    
    substring = strrep(strjoin(regexp(fname{is},'[0-9]','match')),' ','');
    
    FCchans = find(strcmpi({EEG.chanlocs.labels},'FCz'));
    EEG.event(strcmpi({EEG.event.type},'boundary')) = []; % delete boundary events    
    
    % correct for any < 0 SSDs left
    for ie = 1:length(EEG.event)
        if EEG.event(ie).behav(5) < 0; EEG.event(ie).behav(5) = 0; end
        if EEG.event(ie).behav(6) < 0; EEG.event(ie).behav(6) = 0; end
    end
    
    % grab the correct IC comp
    % run algorithmic detection of FC-P3 onset comp, correct sets where the
    % component isn't clear enough for algorithmic detection
    if strcmpi(substring,'10')
        bestcomp = 11;
    elseif strcmpi(substring,'100')
        bestcomp = 3;
    elseif strcmpi(substring,'128') || strcmpi(substring,'233')
        bestcomp = 6;
    elseif strcmpi(substring,'130')
        bestcomp = 13;
    elseif strcmpi(substring,'132') || strcmpi(substring,'16') || strcmpi(substring,'166') || strcmpi(substring,'91') || strcmpi(substring,'77')
        bestcomp = 1;
    elseif strcmpi(substring,'153') || strcmpi(substring,'188') || strcmpi(substring,'224') || strcmpi(substring,'48')
        bestcomp = 5;
    elseif strcmpi(substring,'169') || strcmpi(substring,'192') || strcmpi(substring,'90')
        bestcomp = 2;
    elseif strcmpi(substring,'186')
        bestcomp = 23;
    elseif strcmpi(substring,'189') || strcmpi(substring,'199') || strcmpi(substring,'202') || strcmpi(substring,'22')
        bestcomp = 12;
    elseif strcmpi(substring,'223')
        bestcomp = 10;
    elseif strcmpi(substring,'37')
        bestcomp = 14;
    elseif strcmpi(substring,'54')
        bestcomp = 16;
    elseif strcmpi(substring,'56')
        bestcomp = 4;
    elseif strcmpi(substring,'57')
        bestcomp = 9;
    elseif strcmpi(substring,'96') || strcmpi(substring,'144')
        bestcomp = 7;
    elseif strcmpi(substring,'92')
        bestcomp = 8;
    else
        bestcomp = cipi(EEG);
    end
    % get P3 IC comp
    icweights = 1:size(EEG.icaweights,1);
    
    % extract time series for that comp specifically
    EEG = pop_select(EEG,'nochannel',{'HEOG','VEOG'});    
    EEG = pop_subcomp(EEG, setxor(icweights,bestcomp));
        
    % standardize events so we know which is which
    for ie = 1:length(EEG.event)
        if strcmpi(EEG.event(ie).type,'S  1') % go in go-only
            EEG.event(ie).code = 'Stim';
            if EEG.event(ie).acc == 2  || EEG.event(ie).acc == 99
                EEG.event(ie).type = 'Err/Miss';
            else EEG.event(ie).type = 'CorGo';
            end
        elseif strcmpi(EEG.event(ie).type, 'S  2') % go in stop trial
            EEG.event(ie).code = 'Stim';
            if EEG.event(ie).acc == 3; EEG.event(ie).type = 'SuccStop';
            elseif EEG.event(ie).acc == 4; EEG.event(ie).type = 'FailStop';
            end
        elseif strcmpi(EEG.event(ie).type,'S200') % stop signal
            EEG.event(ie).code = 'StopSignal';
            if EEG.event(ie).acc == 3; EEG.event(ie).type = 'SuccStop';
            elseif EEG.event(ie).acc == 4; EEG.event(ie).type = 'FailStop';
            end
        else EEG.event(ie).type = 'Art'; EEG.event(ie).code = '';
        end
        EEG.event(ie).type = [EEG.event(ie).type '_' EEG.event(ie).code];
    end

    % pool events and associated trial nums
    succstop = EEG.event(~cellfun(@isempty,strfind({EEG.event.type},'SuccStop_StopSignal')));
    succstopgo = EEG.event(~cellfun(@isempty,strfind({EEG.event.type},'SuccStop_Stim')));
    failstop = EEG.event(~cellfun(@isempty,strfind({EEG.event.type},'FailStop_StopSignal')));
    failstopgo = EEG.event(~cellfun(@isempty,strfind({EEG.event.type},'FailStop_Stim')));
    goevents = EEG.event(~cellfun(@isempty,strfind({EEG.event.type},'CorGo_Stim')));
    for ig = 1:length(goevents); goevents(ig).nr = goevents(ig).behav(1); end
    for ig = 1:length(succstop); succstop(ig).nr = succstop(ig).behav(1); end
    for ig = 1:length(failstop); failstop(ig).nr = failstop(ig).behav(1); end
    for ig = 1:length(failstopgo); failstopgo(ig).nr = failstopgo(ig).behav(1); end
    for ig = 1:length(succstopgo); succstopgo(ig).nr = succstopgo(ig).behav(1); end
    % which go signals go with which stop signals
    [a,b,c] = intersect([succstopgo.nr],[succstop.nr]); succstop = succstop(c); succstopgo = succstopgo(b);
    % grab ssds
    for ig = 1:length(succstop)
        if succstop(ig).behav(3) == 1; succstop(ig).ssd = succstop(ig).behav(5); succstopgo(ig).ssd = succstop(ig).ssd; end
        if succstop(ig).behav(3) == 2; succstop(ig).ssd = succstop(ig).behav(6); succstopgo(ig).ssd = succstop(ig).ssd; end
    end
    % same for failed stops
    [a,b,c] = intersect([failstopgo.nr],[failstop.nr]); failstop = failstop(c); failstopgo = failstopgo(b);
    for ig = 1:length(failstop)
        if failstop(ig).behav(3) == 1; failstop(ig).ssd = failstop(ig).behav(5); failstopgo(ig).ssd = failstop(ig).ssd; end
        if failstop(ig).behav(3) == 2; failstop(ig).ssd = failstop(ig).behav(6); failstopgo(ig).ssd = failstop(ig).ssd; end 
    end
    
    trialnumbers = NaN(length(fname),3);
    trialnumbers(is,1) = str2double(is);
    trialnumbers(is,2) = length(succstop);
    trialnumbers(is,3) = length(failstop);
    bldur = round(baseline*EEG.srate/1000);
    epochlength = round(plotlimits*EEG.srate/1000);
    
    % preassign ERP data matrices for gos and stops (with nan, so no chance
    % we're avg'ing in zeros)
    sstopERP = NaN(length(succstop),EEG.nbchan,length(epochlength(1):epochlength(2)));
    fstopERP = NaN(length(failstop),EEG.nbchan,length(epochlength(1):epochlength(2)));
    sgoERP = sstopERP;
    fgoERP = fstopERP;
    
    % for every succ stop
    for ig = 1:length(succstop)
        
        % find the go-only trial that was closest in time after
        nextgo = find([goevents.nr] - succstop(ig).nr > 0,1,'first');
        
        % if there is one, the stop was not the last trial
        if ~isempty(nextgo) && goevents(nextgo).latency + succstop(ig).ssd/2 + epochlength(2) < EEG.pnts
            
            % save all the info about it
            ss_nextgo(ig).num = goevents(nextgo).nr;
            ss_nextgo(ig).behav(8) = goevents(nextgo).behav(8);
            ss_nextgort(ig) = goevents(nextgo).behav(8);
        
            % save the info about the stop
            % ERP: stop signal
            timelock = succstop(ig).latency;
            erpbase = mean(EEG.data(:,timelock+bldur(1):timelock+bldur(2)),2);
            erpbase = repmat(erpbase,1,length(epochlength(1):epochlength(2)));
            sstopERP(ig,:,:) = EEG.data(:,timelock+epochlength(1):timelock+epochlength(2)) - erpbase;
            
            % ERP: matched go
            timelock = goevents(nextgo).latency + succstop(ig).ssd/2;
            erpbase = mean(EEG.data(:,timelock+bldur(1):timelock+bldur(2)),2);
            erpbase = repmat(erpbase,1,length(epochlength(1):epochlength(2)));
            sgoERP(ig,:,:) = EEG.data(:,timelock+epochlength(1):timelock+epochlength(2)) - erpbase;
            
        else
            ss_nextgo(ig).num = NaN;
            ss_nextgo(ig).behav(8) = NaN;
            ss_nextgort(ig) = NaN;  
        end      

    end
    
    % check for and remove any NaN's that popped up as a result of failure
    % to match
    clear remove
    remove = find(isnan(sstopERP(:,1,1)));
    if ~isempty(remove)
        sstopERP(remove,:,:) = [];
        sgoERP(remove,:,:) = [];
        for ir = fliplr(remove')
            succstop(ir) = [];
            succstopgo(ir) = [];
            ss_nextgo(ir) = [];
        end
    end
    
    % do same process for FS
    for ig = 1:length(failstop)

        nextgo = find([goevents.nr] - failstop(ig).nr > 0,1,'first');

        if ~isempty(nextgo) && goevents(nextgo).latency + failstop(ig).ssd/2 + epochlength(2) < EEG.pnts
            fs_nextgo(ig).num = goevents(nextgo).nr;
            fs_nextgo(ig).behav(8) = goevents(nextgo).behav(8);
            fs_nextgort(ig) = goevents(nextgo).behav(8);
        
            % ERP: stop signal
            timelock = failstop(ig).latency;
            erpbase = mean(EEG.data(:,timelock+bldur(1):timelock+bldur(2)),2);
            erpbase = repmat(erpbase,1,length(epochlength(1):epochlength(2)));
            fstopERP(ig,:,:) = EEG.data(:,timelock+epochlength(1):timelock+epochlength(2)) - erpbase;
            
            % ERP: matched go
            timelock = goevents(nextgo).latency + failstop(ig).ssd/2;
            erpbase = mean(EEG.data(:,timelock+bldur(1):timelock+bldur(2)),2);
            erpbase = repmat(erpbase,1,length(epochlength(1):epochlength(2)));
            fgoERP(ig,:,:) = EEG.data(:,timelock+epochlength(1):timelock+epochlength(2)) - erpbase;
        else
            fs_nextgo(ig).num = NaN;
            fs_nextgo(ig).behav(8) = NaN;
            fs_nextgort(ig) = NaN;  
        end

    end
    
    % check for and remove any NaN's that popped up as a result of failure
    % to match
    clear remove
    remove = find(isnan(fstopERP(:,1,1)));
    if ~isempty(remove)
        fstopERP(remove,:,:) = [];
        fgoERP(remove,:,:) = [];
        for ir = fliplr(remove')
            failstop(ir) = [];
            failstopgo(ir) = [];
            fs_nextgo(ir) = [];
        end
    end
       
    % average across trials
    SSTOPERP = squeeze(mean(sstopERP));
    FSTOPERP = squeeze(mean(fstopERP));
    SGOERP = squeeze(mean(sgoERP));
    FGOERP = squeeze(mean(fgoERP));
    SSTOPERR = squeeze(std(sstopERP)/sqrt(length(succstop)));
    FSTOPERR = squeeze(std(fstopERP)/sqrt(length(failstop)));
    SGOERR = squeeze(std(sgoERP)/sqrt(length(succstopgo)));
    FGOERR = squeeze(std(fgoERP)/sqrt(length(failstopgo)));

    % test the SS, FS and their respective matched-go ERPs against each
    % other to produce vectors of p's
    disp(['Running permutations for Subject ' num2str(is) '...']);
    if mcmc == 1
        sample1 = squeeze(mean(sstopERP(:,FCchans,:),2));
        sample2 = squeeze(mean(fstopERP(:,FCchans,:),2));
        sample3 = squeeze(mean(sgoERP(:,FCchans,:),2));
        sample4 = squeeze(mean(fgoERP(:,FCchans,:),2));
        tdiff = zeros(size(sample1,2),1); tdiff_dir = zeros(size(sample1,2),1);
        tdiff2 = zeros(size(sample1,2),1); tdiff_dir2 = zeros(size(sample1,2),1);
        tdiff3 = zeros(size(sample1,2),1); tdiff_dir3 = zeros(size(sample1,2),1);
        for iss = 1:size(sample1,2)
           tdiff2(iss) = mc_ttest(sample3(:,iss),sample1(:,iss),mc);
           tdiff_dir2(iss) = mean(sample1(:,iss)) - mean(sample3(:,iss));
        end
        for iss = 1:size(sample2,2)
           tdiff3(iss) = mc_ttest(sample4(:,iss),sample2(:,iss),mc);
           tdiff_dir3(iss) = mean(sample2(:,iss)) - mean(sample4(:,iss));
        end
    end
    
    % groupdata 
    gERPgsFC(is,:) = SGOERP(FCchans,:);
    gERPsFC(is,:) = SSTOPERP(FCchans,:);
    gERPgfFC(is,:) = FGOERP(FCchans,:);
    gERPfFC(is,:) = FSTOPERP(FCchans,:);
    gERRgsFC(is,:) = SGOERR(FCchans,:);
    gERRsFC(is,:) = SSTOPERR(FCchans,:);
    gERRgfFC(is,:) = FGOERR(FCchans,:);
    gERRfFC(is,:) = FSTOPERR(FCchans,:);
    
    if mcmc == 1; save(fullfile(folders.stats,['Subject' num2str(is) '_P3comp.mat']),'tdiff*'); end
end
% Save group data
save(fullfile(folders.stats,'Group ERP data'),'gERP*','gERR*');

%% Extract P3 onset
p_crit = .05;
SSDlimit = 50; % in ms, trials with SSD below that won't be included
doFDR = 1;
results = zeros(234,4);
visminlat = 1;
searchwindow = [100 200]; % put search window in datapoints (200-400ms)
srate = 500;

for is = 1:length(fname)

    bldur = round(baseline(1)*srate/1000);
    epochlength = round(plotlimits*srate/1000);
    % timelock: select the epoch section that is pre-stim
    timelock = abs(bldur(1)); % no conversion here - epochlength already in datapoints
    
    % load pvals (tdiff: ss>fs; tdiff2: ss>go; tdiff3: fs>go)
    load(fullfile(folders.stats,['Subject' num2str(is) '_P3comp.mat']));    
        
    % This basically finds the peak, then searches backward until the
    % difference is no longer sig. That's onset.
    
    % SUCC: find peak: within search window (excluding the baseline period)
    [~,p3peak] = max(tdiff_dir2(timelock+searchwindow(1):timelock+searchwindow(2)));
    % find peak: within whole epoch (including baseline) in data points
    p3peak = timelock+searchwindow(1)+p3peak;
    % for our own sanity add p3peak in results mat
    results(is,4) = (p3peak-timelock)*1000/srate;
    
    % find ss>go from peak backwards to beginning of epoch
    if doFDR == 1; signif = fdr_bky(tdiff2(timelock:end),p_crit); else; signif = tdiff2(timelock:end)<p_crit; end
    sppp = find(flip(signif(1:p3peak-timelock))==0,1,'first');
    sppp_ms = 1000/srate*(p3peak - sppp - timelock); % FIRST SIGNIFICANT SAMPLE OF SIGN. BLOCK OF SAMPLES TO WHICH PEAK BELONGS
    if ~isempty(sppp_ms); results(is,1) = sppp_ms; else results(is,1) = NaN; end
    
    % find ss>go from peak backwards to beginning of epoch
    % THIS IS THE FIRST DEVIATION FROM STOP TO GO - DOES NOT HAVE TO BE SIG
    pos_pre_peak = find(flip(tdiff_dir2(1:p3peak))<0,1,'first');
    pos_pre_peak = p3peak - pos_pre_peak; % FIRST BEGINNING OF DEFLECTION in samples (from beginning of epoch incl. baseline)
    % put in ms
    ppp_ms = 1000/srate*(pos_pre_peak-timelock); % in ms after stop signal
    if ~isempty(ppp_ms); results(is,2) = ppp_ms; end
    
    % FAIL: find peak: within search window
    [~,p3peak] = max(tdiff_dir3(timelock+searchwindow(1):timelock+searchwindow(2)));
    % find peak: within whole epoch (including baseline) in samples
    p3peak = timelock+searchwindow(1)+p3peak;
    results(is,6) = (p3peak-timelock)*1000/srate;
    % find ss>go from peak backwards to beginning of epoch
    if doFDR == 1; signif = fdr_bky(tdiff3(timelock:end),p_crit); else signif = tdiff3(timelock:end)<p_crit; end
    sppp = find(flip(signif(1:p3peak-timelock))==0,1,'first'); % FIRST SIGNIFICANT SAMPLE OF SIGN. BLOCK OF SAMPLES TO WHICH PEAK BELONGS
    sppp_ms = 1000/srate*(p3peak - sppp - timelock); % FIRST SIGNIFICANT SAMPLE OF SIGN. BLOCK OF SAMPLES TO WHICH PEAK BELONGS
    if ~isempty(sppp_ms); results(is,3) = sppp_ms; else results(is,3) = NaN; end
    
end
%% succ vs. failed ONSET DIFFERENCE ttest
[~,adaptiveonset_p,~,adaptiveonset_stats] = ttest(results(:,1),results(:,3));
adaptiveonset_effect = effectsize(results(:,1),results(:,3));
disp(['P3 onset test: t(' num2str(adaptiveonset_stats.df) ') = ' ...
    num2str(adaptiveonset_stats.tstat) ', p = ' num2str(adaptiveonset_p)...
    ', d = ' num2str(adaptiveonset_effect) '.']);
disp(['On average, SS P3 onsets at ' num2str(mean(results(:,1))) ' and FS P3 at ' num2str(mean(results(:,3))) '.']);

% get max and test peak diff
for is = 1:length(fname)
[SS_AMP(is),SS_PEAK(is)] = max(gERPsFC(is,:));
end
for is = 1:length(fname)
[FS_AMP(is),FS_PEAK(is)] = max(gERPfFC(is,:));
end
[~,peak_p,~,peak_stats] = ttest((SS_PEAK-50)*2,(FS_PEAK-50)*2);
[~,amp_p,~,amp_stats] = ttest(SS_AMP,FS_AMP);
p3_amp_effect = effectsize(SS_AMP,FS_AMP);
p3_peak_effect = effectsize((SS_PEAK-50)*2,(FS_PEAK-50)*2);
disp(['P3 amplitude test: t(' num2str(amp_stats.df) ') = ' ...
    num2str(amp_stats.tstat) ', p = ' num2str(amp_p)...
    ', d = ' num2str(p3_amp_effect) '.']);
disp(['P3 peak timing test: t(' num2str(peak_stats.df) ') = ' ...
    num2str(peak_stats.tstat) ', p = ' num2str(peak_p)...
    ', d = ' num2str(p3_peak_effect) '.']);
disp(['On average, single-sub SS P3 peaks at ' num2str(mean((SS_PEAK-50)*2)) ' and FS P3 at ' num2str((mean(FS_PEAK)-50)*2) '.']);

% N2 amp, peak
for is = 1:length(fname)
N2_SS_AMP(is) = min(gERPsFC(is,125:175)); % look for neg peak b/t 150-250ms
end
for is = 1:length(fname)
N2_FS_AMP(is) = min(gERPfFC(is,125:175)); % look for neg peak b/t 150-250ms
end
[~,N2_amp_p,~,N2_amp_stats] = ttest(N2_SS_AMP,N2_FS_AMP);
n2_amp_effect = effectsize(N2_SS_AMP,N2_FS_AMP);
disp(['N2 amplitude test: t(' num2str(N2_amp_stats.df) ') = ' ...
    num2str(N2_amp_stats.tstat) ', p = ' num2str(N2_amp_p)...
    ', d = ' num2str(n2_amp_effect) '.']);

% P2 amp, peak
for is = 1:length(fname)
P2_SS_AMP(is) = max(gERPsFC(is,100:150)); % look for pos peak b/t 100-200ms
end
for is = 1:length(fname)
P2_FS_AMP(is) = max(gERPfFC(is,100:150)); % look for pos peak b/t 100-200ms
end
[~,P2_amp_p,~,P2_amp_stats] = ttest(P2_SS_AMP,P2_FS_AMP);
p2_amp_effect = effectsize(P2_SS_AMP,P2_FS_AMP);
disp(['P2 amplitude test: t(' num2str(P2_amp_stats.df) ') = ' ...
    num2str(P2_amp_stats.tstat) ', p = ' num2str(P2_amp_p)...
    ', d = ' num2str(p2_amp_effect) '.']);

%% Plot final ERP

% load 
plot(mean(gERPsFC,1),'Color',[30/255,144/255,1],'LineWidth',2); hold on; 
plot(mean(gERPgsFC,1),'Color',[30/255,144/255,1],'LineWidth',2,'LineStyle','--');
plot(mean(gERPfFC,1),'Color',[1,140/255,0],'LineWidth',2); 
plot(mean(gERPgfFC,1),'Color',[1,140/255,0],'LineWidth',2,'LineStyle','--'); 

set(gca,'Ylim',[-2 10]);ax = gca; ax.YLabel.String='Amplitude(mV)'; ax = gca;ax.XLabel.String='Time(ms)'; 
set(gca,'XTick',[0 50 100 150 200 250 300]); set(gca,'XTickLabel',{'-100','SS','100','200','300','400','500'});
set(gca,'Xlim',[0 300]);
line([0 300],[0 0],'Color',[.5 .5 .5]); line([50 50],[-2 10],'Color',[.5 .5 .5]);
legend('Successful stops','Matched go trials (SS)','Failed stops','Matched go trials (FS)');

% add onsets and peaks
line([mean(results(:,1))/2+50 mean(results(:,1))/2+50],[-2 10],'Color',[30/255,144/255,1]);
line([mean(results(:,3))/2+50 mean(results(:,3))/2+50],[-2 10],'Color',[1,140/255,0]);
[~,SS_MAX] = max(mean(gERPsFC,1)); [~,FS_MAX] = max(mean(gERPfFC,1));
line([mean(SS_MAX) mean(SS_MAX)],[-2 10],'Color',[30/255,144/255,1]);
line([mean(FS_MAX) mean(FS_MAX)],[-2 10],'Color',[1,140/255,0]);
