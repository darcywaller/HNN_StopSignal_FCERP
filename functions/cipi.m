% Function for identifying prototypical stop-signal P3 components
function bestcomp = cipi(EEG)

plotlimits = [-.1 .8];
baseline = [-.1 0];

% durations
bldur = baseline*EEG.srate/1000;
epochlength = plotlimits*EEG.srate/1000;

% channels
chans = [find(strcmpi({EEG.chanlocs.labels},'Fz')),find(strcmpi({EEG.chanlocs.labels},'F1')),find(strcmpi({EEG.chanlocs.labels},'F2')),find(strcmpi({EEG.chanlocs.labels},'FCz')),find(strcmpi({EEG.chanlocs.labels},'Cz')),find(strcmpi({EEG.chanlocs.labels},'FC1')),find(strcmpi({EEG.chanlocs.labels},'FC2')),find(strcmpi({EEG.chanlocs.labels},'C1')),find(strcmpi({EEG.chanlocs.labels},'CPz')),find(strcmpi({EEG.chanlocs.labels},'C2')),find(strcmpi({EEG.chanlocs.labels},'CP1')),find(strcmpi({EEG.chanlocs.labels},'CP2'))];
eogchans = [find(strcmpi({EEG.chanlocs.labels},'VEOG')),find(strcmpi({EEG.chanlocs.labels},'HEOG'))];
EEG = pop_select(EEG,'nochannel',eogchans);

% topographical criterion
icawinv = EEG.icawinv; % the ica weights
tc = 0; top = [];
for ic = 1:size(icawinv,2)
    %[~,topo] = outlier_z(icawinv(:,ic),pcrit);
    [~,topo(1)] = max(icawinv(:,ic));% check for mins and maxs at our selected electrodes
    [~,topo(2)] = min(icawinv(:,ic));
    % hold as candidate if the max activation at any specified chans
    if ~isempty(intersect(topo,chans)); tc = tc+1; top(tc) = ic; end
end

% load comperps
csERP = zeros(length(ic), 450);

for ic = 1:length(top) % cycle through the comps that met topo filters
    icweights = [1:size(EEG.icaweights,1)];
    cEEG = pop_subcomp(EEG, setxor(icweights,top(ic))); % get time series for comp
    cEEG = pop_epoch(cEEG,{'S200'},plotlimits,'epochinfo', 'yes');
    cEEG = pop_rmbase(cEEG,[-100,0]);
    % make ERP for SS, FS
    stopEEG = pop_selectevent(cEEG,'type',{'S200'},'acc',3,'deleteevents','off','deleteepochs','on','invertepochs','off');
    failEEG = pop_selectevent(cEEG,'type',{'S200'},'acc',4,'deleteevents','off','deleteepochs','on','invertepochs','off');
    
    erp.settings.electrode = [find(strcmpi({EEG.chanlocs.labels},'Fz')),find(strcmpi({EEG.chanlocs.labels},'FCz')),find(strcmpi({EEG.chanlocs.labels},'Cz')),find(strcmpi({EEG.chanlocs.labels},'C1')),find(strcmpi({EEG.chanlocs.labels},'C2'))];
    electrodedata_stop = stopEEG.data(erp.settings.electrode,:,:);
    meanelectrodedata_stop = nanmean(electrodedata_stop,3);
    electrodedata_fail = failEEG.data(erp.settings.electrode,:,:);
    meanelectrodedata_fail = nanmean(electrodedata_fail,3);
    % take mean diff between SS and FS ERP (which will involve onset diffs)
    csERP(ic,:) = nanmean(meanelectrodedata_stop,1)-nanmean(meanelectrodedata_fail,1);
    
end


% make stop erps for this subject overall
EEG = pop_epoch(EEG,{'S200'},plotlimits,'newname','Merged datasets pruned with ICA epochs', 'epochinfo', 'yes');
EEG = pop_rmbase(EEG, [-50,0]);
stopEEG = pop_selectevent(EEG,'type',{'S200'},'acc',3,'deleteevents','off','deleteepochs','on','invertepochs','off');
electrodedata_stop = stopEEG.data(erp.settings.electrode,:,:);
sstoperp = mean(nanmean(electrodedata_stop,3),1);
failEEG = pop_selectevent(EEG,'type',{'S200'},'acc',4,'deleteevents','off','deleteepochs','on','invertepochs','off');
electrodedata_fail = failEEG.data(erp.settings.electrode,:,:);
fstoperp = mean(nanmean(electrodedata_fail,3),1);
diffstoperp = zeros(1,length(sstoperp));
% take mean diff between SS and FS ERP (which will involve onset diffs)
for ip = 1:length(sstoperp)
    diffstoperp(ip) = sstoperp(ip) - fstoperp(ip);
end

% correlations between comp ERPs and all-IC ERPs
window = [150 300]; % 200-500ms
R = zeros(length(top), 3);
for ic = 1:size(csERP,1)
    [r,p] = corrcoef(diffstoperp(1,window(1):window(2)),csERP(ic,window(1):window(2)));
    % save correlation coeff, p, and which comp
    R(ic,1) = r(1,2); R(ic,2) = p(1,2); R(ic,3) = top(ic);
end
[~,bestcor] = max(R(:,1));

bestcomp = R(bestcor,3);
