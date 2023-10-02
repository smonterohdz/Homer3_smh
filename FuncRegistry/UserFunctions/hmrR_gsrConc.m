% SYNTAX:
% [dc] = hmrR_gsrConc(dc, probe, mlActMan, mlActAuto, rhoSD_ssThresh)
%
% UI NAME:
% hmrR_gsrConc
%
% DESCRIPTION:
% This script
% Global signal regression 
%
% INPUT:
% dc - SNIRF data type where dataTimeCourse is concentration data (See SNIRF Spec for more details)
% probe - SNIRF probe type containing source/detector geometry data (See SNIRF Spec for more details)
% mlActMan: Cell array of vectors, one for each time base in data, specifying
%        active/inactive channels with 1 meaning active, 0 meaning inactive
% mlActAuto: Cell array of vectors, one for each time base in data, specifying
%        active/inactive channels with 1 meaning active, 0 meaning inactive
% rhoSD_ssThresh - max distance for a short separation measurement. 
%
%
% OUTPUTS:
% dc - SNIRF.data container with the Hb with global mean subtracted.
%
% USAGE OPTIONS:
% hmrR_gsrConc: [dc] = hmrR_gsrConc(dc, probe, mlActMan, mlActAuto, rhoSD_ssThresh)
%
%
% PARAMETERS:
% rhoSD_ssThresh: 15.0
%
% PREREQUISITES:
% For the hmrR_gsrConc Usage Option, use Delta_OD_to_Conc: dc = hmrR_OD2Conc( dod, probe, ppf )
%
%
function [dc] = hmrR_gsrConc(dc, probe, mlActMan, mlActAuto, rhoSD_ssThresh)

% initialize active channel list
if isempty(mlActMan)
    mlActMan = cell(length(dc),1);
end
if isempty(mlActAuto)
    mlActAuto = cell(length(dc),1);
end

SrcPos  = probe.GetSrcPos();
DetPos  = probe.GetDetPos();

for iBlk=1:length(dc)

    % Extract variables from array-of-data-blocks SNIRF arguments, hence the processing
    % from within the for loop. (see SNIRF spec on github.com for details on data blocks)
    d        = dc(iBlk).GetDataTimeSeries('reshape');
    datatype = dc(iBlk).GetDataTypeLabel();  % Get the input data type
    t        = dc(iBlk).GetTime();
    ml       = dc(iBlk).GetMeasListSrcDetPairs('reshape');

    fq = 1/(t(2)-t(1));

    % get  a list of active channels
    mlActMan{iBlk} = mlAct_Initialize(mlActMan{iBlk}, ml);
    mlActAuto{iBlk} = mlAct_Initialize(mlActAuto{iBlk}, ml);
    MeasListAct = mlActMan{iBlk}(:,3) & mlActAuto{iBlk}(:,3);

    %% find the list of short and long distance channels
    lst = 1:size(ml,1);
    rhoSD = zeros(length(lst),1);
    posM = zeros(length(lst),3);
    for iML = 1:length(lst)
        rhoSD(iML) = sum((SrcPos(ml(lst(iML),1),:) - DetPos(ml(lst(iML),2),:)).^2).^0.5;
        posM(iML,:) = (SrcPos(ml(lst(iML),1),:) + DetPos(ml(lst(iML),2),:)) / 2;
    end
    lstSS = lst(find(rhoSD<=rhoSD_ssThresh &  MeasListAct(lst)==1)); %#ok<*FNDSB>
    lstLS = lst(find(rhoSD>rhoSD_ssThresh & MeasListAct(lst)==1));

    %% get long and short separation data
    if strncmp(datatype{1}, 'Hb', 2)
        dHbO = squeeze(d(:,1,:));
        dHbR = squeeze(d(:,2,:));
        d_short = [dHbO(:,lstSS), dHbR(:,lstSS)];
        d_long  = [dHbO(:,lstLS), dHbR(:,lstLS)];
    else
        return;
    end

    %% GSR
    y = dHbO;
    yMean = mean(y(:,lstSS),2);
    g = [ones(size(yMean)) yMean];
    % Regularization parameter
    lambda = 0.1*max(diag(g'*g));
    % Regressor: Tikhonov regression estimator
    betaCoeff = pinv(g'*g + lambda*eye(size(g'*g)))*g'*y;
    % Filtered signal yHat
    d(:,1,:) = y - g*betaCoeff;

    y = dHbR;
    yMean = mean(y(:,lstSS),2);
    g = [ones(size(yMean)) yMean];
    % Regularization parameter
    lambda = 0.1*max(diag(g'*g));
    % Regressor: Tikhonov regression estimator
    betaCoeff = pinv(g'*g + lambda*eye(size(g'*g)))*g'*y;
    % Filtered signal yHat
    d(:,2,:) = y - g*betaCoeff;
    
    %%
    dc.SetDataTimeSeries(d);

    %     s = SnirfClass();
    %     s.data = dc;
    %     s.aux = aux;
    %     AUX = s.GetAuxDataMatrix();
    %     AUX = AUX(:,tCCAaux_inx);
end