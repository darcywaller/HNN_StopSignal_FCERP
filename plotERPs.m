%% Plot ERPs and report stats from pre-extracted data
% Written by Darcy Diesburg (2023), Brown University

% some settings
baseline = [-100 0];
plotlimits = [-100 1000];
%% Quantify P3 onset at the single-subject level
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
load(fullfile(pwd,'/stats/Group ERP data.mat'));
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