%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSRT >> SSRT computation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUT (all inputs have to be n*1 vectors of the same size and order):
%           
%       RT          = all RTs (Go, failed inhibitit, succesfull inhibit)
%       acc         = Outcomes for stoptrials (1 = success, 2 = failure 0 = go)
%       SSD         = Stop signal delays for stoptrials
%   
%   OUTPUT
%
%       SSRT_mean   = SSRT based on the standard mean method
%       SSRT_cmean  = SSRT based on the convergence method
%                        (criterion by I. Greenhouse)
%       convergence = points of convergence for every staircase
%       p_inhibit   = p(inh) for every staircase from its point of convergence
%       SSRT_integ  = SSRT based o the integration method
%
% Jan R. Wessel, University of California, San Diego, October 2011
%   Email: jwessel@ucsd.edu
%
%   programmed on MATLAB 2009a
%
%   Version History:
%       11/22/2011: Version 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SSRTdata = SSRT(RT, acc, SSD, staircase)

% check inputs: nargin
if nargin > 5
    error('SSRT.m: Too many input variables. Type "help SSRT" for more info');
elseif nargin < 1
    help SSRT
    return
elseif nargin < 4
    disp('SSRT.m: Using default value for staircase: 1');
    staircase = ones(length(SSD),1);
elseif nargin < 3
    error('SSRT.m: Too few input variables. Provide at least RT, acc, and SSD. Type "help SSRT" for more info');
end
% check inputs: dimensions
if min(size(RT)) > 1 || min(size(acc)) > 1 || min(size(SSD)) > 1
    error('SSRT.m: One of the variables is not a vector.');
end
if size(RT,1) ~= size(acc,1) || size(RT,1) ~= size(SSD,1) || size(acc,1) ~= size(SSD,1) 
    error('SSRT.m: The input vectors have different dimensionality.');
end

% combine inputs
matrix = [RT acc SSD staircase];

%% mean method
% SSD mean
scs = unique(staircase(staircase>0)); % get all staircases
for is = 1:numel(scs)
    % get SSD mean of all staircases
    sctrials = matrix(matrix(:,4)==is,:);
    SSRTdata.scSSD(is) = nanmean(sctrials(:,3));
end
SSRTdata.meanSSD = mean(SSRTdata.scSSD);
% go trial mean
gos = matrix(matrix(:,2)==0,:); % get all trials for acc = 0, which are the go trials
SSRTdata.meanGoRT = nanmean(gos(:,1));
% average
SSRTdata.SSRT_mean = SSRTdata.meanGoRT - SSRTdata.meanSSD;

%% convergence method
% for each staircase, get point of convergence based on ians criterion
for is = 1:numel(scs)
    % get trials for each staircase
    sctrials = matrix(matrix(:,4)==is,:);
    if size(sctrials,1) < 5 % too few trials in staircase
        disp('SSRT.m Too few trials in analysis,setting NaN');
        SSRTdata.scSSD_conv(is) = NaN;
        break
    end
    % get +/- ms
    SSRTdata.difference(is) = abs(sctrials(1,3) - sctrials(2,3));
    % find convergence according to +-+ sequence (Ian)
    differentiation = diff(sctrials(:,3));
    for it = 1:length(differentiation)-2 % integrate over three trials and find point where this is +difference
        nextseq(it) = sum(differentiation(it:it+2));
    end
    if ~isempty(find(nextseq == SSRTdata.difference(is),1,'first'))
        SSRTdata.convergence(is) = find(nextseq == SSRTdata.difference(is),1,'first');
        
        % Get trials after convergence
        convtrials = sctrials(SSRTdata.convergence(is):end,:);
        % get p_inhibit
        SSRTdata.p_inhibit(is) = sum(convtrials(:,2)==1) / (sum(convtrials(:,2)==2) + sum(convtrials(:,2)==1));

        % get mean SSD from convergence point to end for each staircase
        SSRTdata.scSSD_conv(is) = nanmean(convtrials(:,3));
    else
        disp('No convergence, skipping convergence part and setting NaN');
        SSRTdata.scSSD_conv(is) = NaN;
    end
    
end
% get SSRT
SSRTdata.SSRT_conv = SSRTdata.meanGoRT - mean(SSRTdata.scSSD_conv);

%% Integration method
SSRT_int = nan(numel(scs),1);
for isc = 1:numel(scs) % for each staircase
    
    % get unique SSDs
    uSSD = unique(matrix(:,3));
    % get number of trials and p(inh) at each value
    for is = 1:numel(uSSD)
        SSDtrials = matrix(matrix(:,3) == uSSD(is),:);
        n_SSD(is) = size(SSDtrials,1);
        p_inh(is) = sum(SSDtrials(:,2)==1) / (sum(SSDtrials(:,2)==2) + sum(SSDtrials(:,2)==1));
    end
    % make matrix to sort and go thru
    SSDmatrix = [uSSD n_SSD' p_inh'];
    SSDmatrix = flipud(sortrows(SSDmatrix,[2 3])); % sort by number of trials, then p inhibit, then flip
    % delete zero prob lines
    SSDmatrix(SSDmatrix(:,3)<=0,:) = [];
    SSDmatrix(isnan(SSDmatrix(:,3)),:) = [];
    
    % get individual SSRTs at most frequent SSDs
    SSRT_per_SSD = zeros(size(SSDmatrix,1),1);
    for is = 1:size(SSDmatrix,1)
        SSRT_per_SSD(is) = quantile(gos(:,1),(1 - SSDmatrix(is,3))) - SSDmatrix(is,1); % get p(inh) quantile in go RT dist and subtract SSD value
    end
    
    if length(SSRT_per_SSD)>1
    % integrate across most frequent 2 SSDs
        SSRT_int(isc) = mean(SSRT_per_SSD(1:2));
    else
        SSRT_int(isc) = SSRT_per_SSD(1);
    end
    
end
% average across staircases
SSRTdata.SSRT_integ = mean(SSRT_int);