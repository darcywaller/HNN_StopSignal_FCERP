function behav = Basicstop_analyze(trialseq,replace_misses)

if nargin < 2; replace_misses = 1; end

id = Basicstop_columns;

% see if there are subzero ssds
trialseq(trialseq(:,id.lssd) < 0,id.lssd) = 0; 
trialseq(trialseq(:,id.rssd) < 0,id.rssd) = 0; 

% trialgroups
gotrials = trialseq(trialseq(:,id.stop)==0,:);
stoptrials = trialseq(trialseq(:,id.stop)==1,:);
corgotrials = gotrials(gotrials(:,id.acc)==1,:);
corgotrialsL = corgotrials(corgotrials(:,id.dir)==1,:);
corgotrialsR = corgotrials(corgotrials(:,id.dir)==2,:);
errtrials = gotrials(gotrials(:,id.acc)==2,:);
succstoptrials = stoptrials(stoptrials(:,id.acc)==3,:);
failstoptrials = stoptrials(stoptrials(:,id.acc)==4,:);
misstrials = gotrials(gotrials(:,id.resp)==0,:);
% for ssrt (this is an edit makde in January 2019 to implement the
% go-failure replacement method)
ssrttrials = trialseq(trialseq(:,id.acc)~=2,:); % all but direction errors
ssrttrials(ssrttrials(:,id.acc) == 5 | ssrttrials(:,id.acc) == 99,id.rt) = max(corgotrials(:,id.rt)); % decided below if they will be included or not

% post-trials
[~,postss] = intersect(corgotrials(:,id.num),succstoptrials(:,id.num)+1);
[~,postfs] = intersect(corgotrials(:,id.num),failstoptrials(:,id.num)+1);
[~,postgo] = intersect(corgotrials(:,id.num),corgotrials(:,id.num)+1);
postsstrials = corgotrials(postss,:);
postfstrials = corgotrials(postfs,:);
postgotrials = corgotrials(postgo,:);

% trialnumbers
behav.numbers.go.all = size(gotrials,1);
behav.numbers.go.cor = size(corgotrials,1);
behav.numbers.go.err = size(errtrials,1);
behav.numbers.go.miss = size(misstrials,1);
behav.numbers.stop.all = size(stoptrials,1);
behav.numbers.stop.succ = size(succstoptrials,1);
behav.numbers.stop.fail = size(failstoptrials,1);

% rates
behav.rates.error = behav.numbers.go.err / behav.numbers.go.all;
behav.rates.miss = behav.numbers.go.miss / behav.numbers.go.all;
behav.rates.stop = behav.numbers.stop.all / size(trialseq,1);
behav.rates.stopsucc = behav.numbers.stop.succ / behav.numbers.stop.all;

% RT
behav.RT.corgo = mean(corgotrials(:,id.rt));
behav.RT.corgoL = mean(corgotrialsL(:,id.rt));
behav.RT.corgoR = mean(corgotrialsR(:,id.rt));
behav.RT.corgo_median = median(corgotrials(:,id.rt));
behav.RT.errgo = mean(errtrials(:,id.rt));
behav.RT.failstop = mean(failstoptrials(:,id.rt));
behav.RT.postgo = mean(postgotrials(:,id.rt));
behav.RT.postss = mean(postsstrials(:,id.rt));
behav.RT.postfs = mean(postfstrials(:,id.rt));

% SSRT
ssrtmat = zeros(size(ssrttrials,1),4);
ssrtmat = [ssrttrials(:,id.rt) ssrttrials(:,id.acc)];
for it = 1:size(ssrttrials,1)
    ssrtmat(it,3) = ssrttrials(it,id.stop+ssrttrials(it,id.dir)); %SSD
    ssrtmat(it,4) = ssrttrials(it,id.dir);
    ssrtmat(it,5) = ssrttrials(it,id.blo);
end
ssrtmat(ssrtmat(:,3) == 0,3) = 0;
if replace_misses == 1 % misses will be counted towards gort (with rt = max(rt))
    ssrtmat(ssrtmat(:,2) == 1 | ssrtmat(:,2) == 5 | ssrtmat(:,2) == 99,2) = 0;
else; ssrtmat(ssrtmat(:,2) == 1,2) = 0; % misses are still in ssrtmat as 5/99 but will not be counted by SSRT function
end
ssrtmat(ssrtmat(:,2) == 3,2) = 1;
ssrtmat(ssrtmat(:,2) == 4,2) = 2;
ssrtmat(ssrtmat(:,2) == 0,4) = 0;
SSRTdata = SSRT(ssrtmat(:,1), ssrtmat(:,2), ssrtmat(:,3), ssrtmat(:,4));
behav.RT.SSRT = SSRTdata.SSRT_mean;
behav.RT.SSRTi = SSRTdata.SSRT_integ;
behav.RT.mSSD = SSRTdata.meanSSD;
