% SYNTAX:
% mlActAuto = hmrR_PruneChannelsQTNIRS(data, probe, mlActMan, tIncMan, windowSec,hpf,lpf,ScanQualityThresh,SCIthresh,PSPthresh)

%
% UI NAME:
% PruneChannelsQTNIRS
%
% DESCRIPTION:
% Prune channels from the measurement list based on the presence of cardiac 
% pulsation. Using the combination of SCI and PSP metrics (and their
% thresholds), this function updates MeasListAct based on whether data 'd' 
% meets these conditions as specified by 'SCIthreshold', 'PSPthreshold', 
% and 'QualityThreshold'.
%
% INPUTS:
% data - SNIRF object containing time course data (nTpts x nChannels )
% probe - SNIRF object describing the probe - optode positions and wavelengths.
% mlActMan - manually selected active channels
% windowSec - Length (sec.) of the windows to compute SCI and PSP
% hpf - Max frequency for the band-pass filtering
% lpf - Min frequency for the band-pass filtering
% ScanQualityThresh - The required quality value (normalized; 0-1) of good-quality windows in every channel.
% SCIthresh - Minimum value accepted for SCI (0-1)
% PSPthresh - Minimum value for PSP (0-0.1)

%
% OUTPUTS:
% mlAct - cell array of all data blocks - each data block is an array
%         of true/false for all channels in the contanining data block
%         specifying active/inactive status. (# of data blocks x # of Channels)
%
% USAGE OPTIONS:
% PruneChannelsQTNIRS: mlActAuto = hmrR_PruneChannelsQTNIRS(data, probe, mlActMan, tIncMan, windowSec,hpf,lpf,ScanQualityThresh,SCIthresh,PSPthresh)
%
% PARAMETERS:
% windowSec: 5
% hpf: 0.5
% lpf: 2.5
% ScanQualityThresh: 0.7
% SCIThresh: 0.7
% PSPThresh: 0.08
%
% PREREQUISITES:
% For the hmrR_PruneChannelsQTNIRS Usage Option, use Intensity values.
%
% TO DO:
% 
%
function mlActAuto = hmrR_PruneChannelsQTNIRS(data, probe, mlActMan, tIncMan, windowSec,hpf,lpf,ScanQualityThresh,SCIthresh,PSPthresh)

% Init output 
mlActAuto = cell(length(data),1);

% Check input args
if nargin<10
    disp( 'USAGE: hmrR_PruneChannelsQTNIRS(data, probe, mlActMan, tIncMan, windowSec,hpf,lpf,ScanQualityThresh,SCIthresh,PSPthresh)' )
    return
end
if isempty(tIncMan)
    tIncMan = cell(length(data),1);
end
if isempty(mlActMan)
    mlActMan = cell(length(data),1);
end

for iBlk = 1:length(data)

    raw.d        = data(iBlk).GetDataTimeSeries();
    raw.t        = data(iBlk).GetTime();
    raw.SD.MeasList     = data(iBlk).GetMeasList();
    raw.SD.Lambda       = probe.GetWls();
    raw.SD.SrcPos = probe.sourcePos3D;
    raw.SD.DetPos = probe.detectorPos3D;
   
    raw.s        = zeros(size(raw.t));

    SrcPos   = probe.GetSrcPos();
    DetPos   = probe.GetDetPos();

    mlActMan{iBlk} = mlAct_Initialize(mlActMan{iBlk}, raw.SD.MeasList);

    if isempty(tIncMan{iBlk})
        tIncMan{iBlk} = ones(length(raw.t),1);
    end
    %tInc = tIncMan{iBlk};

    %lstInc = find(tInc==1);
    %d = d(lstInc,:);

    fcut = [hpf lpf];
    quality_threshold = ScanQualityThresh;
    sci_threshold = SCIthresh;
    psp_threshold = PSPthresh;

    %lstInc = find(tInc==1);
    %d = d(lstInc,:);
 

    %Sortting rows by source, detector, and lambda
    %[raw.SD.MeasList, idxML] = sortrows(raw.SD.MeasList,[1 2 4]);
    %raw.d = raw.d(:,idxML);


    if ~exist('fcut','var')
        fcut = [0.5 2.5];
    end
    if ~exist('windowSec','var')
        windowSec = 5;
    end
    if ~exist('windowOverlap','var')
        windowOverlap = 0;
    end
    if ~exist('quality_threshold','var')
        quality_threshold = 0.5;
    end
    if ~exist('sci_threshold','var')
        sci_threshold = 0.7;
    end
    if ~exist('psp_threshold','var')
        psp_threshold = 0.08;
    end
    if ~exist('cond_mask','var')
        cond_mask = 'resting';
    end
    if ~exist('lambda_mask_','var')
        lambdas_ = unique(raw.SD.MeasList(:,4))';
        lambdas_ = lambdas_(1:2);
        lambda_mask = true(1,length(lambdas_));
        if length(lambda_mask) ~= length(raw.SD.Lambda)
            for ii=3:length(raw.SD.Lambda)
                lambda_mask(end+1) = 0;
            end
        end
    end
    if ~exist('dodFlag','var')
        dodFlag = -1;
    end
    if ~exist('guiFlag','var')
        guiFlag = 0;
    end

    qualityMats = qtnirs(raw,'freqCut',fcut,...
        'window',windowSec,...
        'overlap',windowOverlap,....
        'qualityThreshold',quality_threshold,...
        'sciThreshold',sci_threshold,...
        'pspThreshold',psp_threshold,...
        'conditionsMask',cond_mask,...
        'dodFlag',dodFlag,...
        'guiFlag',guiFlag);
    % update SD.MeasListAct

    % Start by including all channels
    chanList = ones(size(raw.SD.MeasList,1),1);
    
    idxqmWL1 = qualityMats.MeasList(:,4)==1;
    idxSDWL1 = raw.SD.MeasList(:,4)==1;
    reswl1=qualityMats.MeasListAct(idxqmWL1);
    chanList(idxSDWL1) = reswl1;

    idxqmWL2 = qualityMats.MeasList(:,4)==2;
    idxSDWL2 = raw.SD.MeasList(:,4)==2;
    reswl2=qualityMats.MeasListAct(idxqmWL2);
    chanList(idxSDWL2) = reswl2;
    
    chanList = chanList(:) & mlActMan{iBlk}(:,3);

    % update MeasListAct  
    mlActAuto{iBlk} = mlAct_Initialize(chanList, raw.SD.MeasList);


    % --------------------------------------------------------------------
    % 
    % 
    % %--------------------------------
    % % check for dRange and SNRthresh
    % dmean = mean(d,1);
    % dstd = std(d,[],1);
    % 
    % idxs1 = [];
    % idxs2 = [];
    % idxs3 = [];
    % nLambda = length(Lambda);
    % lst1 = find(MeasList(:,4)==1);    
    % 
    % % Start by including all channels
    % chanList = ones(size(MeasList,1),1);
    % 
    % for ii = 1:nLambda
    %     lst = [];
    %     rhoSD = [];
    %     for jj = 1:length(lst1)
    %         lst(jj) = find(MeasList(:,1)==MeasList(lst1(jj),1) & ...
    %                        MeasList(:,2)==MeasList(lst1(jj),2) & ...
    %                        MeasList(:,4)==ii);
    %         rhoSD(jj) = norm( SrcPos(MeasList(lst1(jj),1),:) - DetPos(MeasList(lst1(jj),2),:) );
    %     end
    % 
    %     % dRange exclusion criteria
    %     idxs1 = [idxs1, find(dmean(lst)<=dRange(1) | dmean(lst)>=dRange(2))];
    % 
    %     % SNRthresh exclusion criteria
    %     idxs2 = [idxs2, find((dmean(lst)./dstd(lst)) <= SNRthresh)];
    % 
    %     % SDrange exclusion criteria
    %     idxs3 = [idxs3, find(rhoSD<SDrange(1) | rhoSD>SDrange(2))];
    % 
    %     idxsExcl = unique([idxs1(:)', idxs2(:)', idxs3(:)']);
    % 
    %     chanList(lst(idxsExcl)) = 0;
    % end
    % 
    % idxs4 = find(mlActMan{iBlk}(:,3) == 0);
    % 
    % if ~isempty(idxs1)
    %     fprintf('hmrR_PruneChannels:  excluded channels [ %s ] based on dRange=[ %s ] at wavelength %d\n', num2str(idxs1(:)'), num2str(dRange), ii);
    % end
    % if ~isempty(idxs2)
    %     fprintf('hmrR_PruneChannels:  excluded channels [ %s ] based on SNRthresh=%0.1f at wavelength %d\n', num2str(idxs2(:)'), SNRthresh, ii);
    % end
    % if ~isempty(idxs3)
    %     fprintf('hmrR_PruneChannels:  excluded channels [ %s ] based on SDrange=[ %s ] at wavelength %d\n', num2str(idxs3(:)'), num2str(SDrange), ii);
    % end
    % if ~isempty(idxs4)
    %     fprintf('hmrR_PruneChannels:  manually excluded channels [ %s ]\n', num2str(idxs4(:)'));
    % end
    % 
    % chanList = chanList(:) & mlActMan{iBlk}(:,3);
    % %---------------------------------------------------
    % 
    % % update MeasListAct  
    % mlActAuto{iBlk} = mlAct_Initialize(chanList, MeasList);
end

