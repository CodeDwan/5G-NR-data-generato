function [] = datagenerator_rayt(bsPosition,bsHight,bsrange,mapf,savefile,samplenum)


format long
bsset.bsPosition = bsPosition; % Lat, lon
bsset.bsHight =bsHight;                       
bsset.bsrange = bsrange;   
mapfile = mapf;       %选择地图

filename = savefile;
%% 参数设置
number = samplenum; 
freq = 28e9;  
RBNum = 26;   %资源块数
BSTX = [8,4];   %基站天线
UETX = [2,2];    %用户天线
SCS = 30; % subcarrier spacing
NRB = RBNum; % number of resource blocks, 10 MHz bandwidth
HCSI = zeros(samplenum,RBNum, prod(UETX),prod(BSTX));
%% 其他参数
bsset.bsAntSize = BSTX;                    
bsset.bsArrayOrientation = [0 0].';       % 
degrange = bsset.bsrange / 111000; 

userset.lat = randpostion(bsset.bsPosition(1)-degrange,bsset.bsPosition(1)+degrange,number);
userset.lont = randpostion(bsset.bsPosition(2)-degrange,bsset.bsPosition(2)+degrange,number);
userset.hight = randpostion(0.5,1.5,number);
userset.azi = randpostion(0,360,number);
userset.ele = randpostion(0,90,number);
userset.Ant = UETX;

for index = 1:number
    % 该样本参数设置_bs
    fc = freq;                            % carrier frequency (Hz)
    bsPosition = bsset.bsPosition;
    bsHight = bsset.bsHight;
    bsArrayOrientation = bsset.bsArrayOrientation;
    bsAntSize = bsset.bsAntSize;
    % 该样本参数设置_ue
    uePosition = [userset.lat(index),userset.lont(index)];  % Lat, lon
    uehight = userset.hight(index);
    ueAntSize = userset.Ant;                    % number of rows and columns in rectangular array (UE).
    ueArrayOrientation = [userset.azi(index) userset.ele(index)].';      % azimuth (0 deg is East, 90 deg is North) and elevation (positive points upwards)  in deg
    reflectionsOrder = 1;                 % number of reflections for ray tracing analysis (0 for LOS)

    %% Import and Visualize 3-D Environment with Buildings for Ray Tracing
    % Launch Site Viewer with buildings in Hong Kong. For more information about 
    % the osm file, see [1].

    if exist('viewer','var') && isvalid(viewer) % viewer handle exists and viewer window is open
        viewer.clearMap();
    else
        viewer = siteviewer("Basemap","openstreetmap","Buildings",mapfile); 
        % viewer = siteviewer("Basemap","darkwater","Buildings","hongkong.osm"); 
    end

    %% 
    % 
    %% Create Base Station and UE
    % Locate the base station and the UE on the map.

    bsSite = txsite("Name","Base station", ...
        "Latitude",bsPosition(1),"Longitude",bsPosition(2),...
        "AntennaAngle",bsArrayOrientation(1:2),...
        "AntennaHeight",bsHight,...  % in m
        "TransmitterFrequency",fc);

    ueSite = rxsite("Name","UE", ...
        "Latitude",uePosition(1),"Longitude",uePosition(2),...
        "AntennaHeight",uehight,... % in m
        "AntennaAngle",ueArrayOrientation(1:2));

    
    bsSite.show();
    ueSite.show();
   
    %% Ray Tracing Analysis
    % Perform ray tracing analysis using the method of images. This method models 
    % surface specular reflections but does not include effects from refraction, diffraction, 
    % scattering, or transmission through buildings. 

    rays = raytrace(bsSite,ueSite,"NumReflections",0:reflectionsOrder,"Type","pathloss");
    
    if isempty(rays{1})
        continue
    else
        plot(rays{1})

        pathToAs = [rays{1}.PropagationDelay]-min([rays{1}.PropagationDelay]);  % Time of arrival of each ray (normalized to 0 sec)
        avgPathGains  = -[rays{1}.PathLoss];                                    % Average path gains of each ray
        pathAoDs = [rays{1}.AngleOfDeparture];                                  % AoD of each ray
        pathAoAs = [rays{1}.AngleOfArrival];                                    % AoA of each ray
        isLOS = any([rays{1}.LineOfSight]);                                     % Line of sight flag

        channel = nrCDLChannel;
        channel.DelayProfile = 'Custom';
        channel.PathDelays = pathToAs;
        channel.AveragePathGains = avgPathGains;
        channel.AnglesAoD = pathAoDs(1,:);       % azimuth of departure
        channel.AnglesZoD = 90-pathAoDs(2,:);    % channel uses zenith angle, rays use elevation
        channel.AnglesAoA = pathAoAs(1,:);       % azimuth of arrival
        channel.AnglesZoA = 90-pathAoAs(2,:);    % channel uses zenith angle, rays use elevation
        channel.HasLOSCluster = isLOS;
        channel.CarrierFrequency = fc;
        channel.NormalizeChannelOutputs = false; % do not normalize by the number of receive antennas, this would change the receive power
        channel.NormalizePathGains = false;      % set to false to retain the path gains

        c = physconst('LightSpeed');
        lambda = c/fc;

        % UE array
        ueArray = phased.URA('Size',ueAntSize(1:2),'ElementSpacing', 0.5*lambda*[1 1]);            % isotropic element by default
        channel.ReceiveAntennaArray = ueArray;
        channel.ReceiveArrayOrientation = [ueArrayOrientation(1); (-1)*ueArrayOrientation(2); 0];  % the (-1) converts elevation to downtilt

        % Base station array
        bsArray = phased.URA('Size',bsAntSize(1:2),'ElementSpacing', 0.5*lambda*[1 1],'Element',phased.NRAntennaElement);
        channel.TransmitAntennaArray = bsArray;
        channel.TransmitArrayOrientation = [bsArrayOrientation(1); (-1)*bsArrayOrientation(2); 0];   % the (-1) converts elevation to downtilt


        ofdmInfo = nrOFDMInfo(NRB,SCS);

        channel.SampleRate = ofdmInfo.SampleRate;
        %% Channel Estimation

        channel.ChannelFiltering = false;
        [pathGains,sampleTimes] = channel();


        pg=permute(pathGains,[2 1 3 4]); % first dimension is the number of paths
        if isLOS
            % in LOS cases sum the first to paths, they correspond to the LOS ray
            pg = [sum(pg(1:2,:,:,:)); pg(3:end,:,:,:)];
        end
        pg = abs(pg).^2;

        %% 
        % Obtain a perfect channel estimate for slot 0.

        pathFilters = getPathFilters(channel);
        nSlot = 0;
        [offset,~] = nrPerfectTimingEstimate(pathGains,pathFilters);
        if offset > 35
            continue
        else
            hest = nrPerfectChannelEstimate(pathGains,pathFilters,NRB,SCS,nSlot,offset,sampleTimes);
        end

        %% 
        HCSI(index,:,:,:) = dealCSI(hest,RBNum);
        %%

        surf(pow2db(abs(hest(:,:,1,1)).^2));
        shading('flat');
        xlabel('OFDM Symbols');ylabel('Subcarriers');zlabel('Magnitude Squared (dB)');
        title('Channel Magnitude Response (1^{st} tx - 1^{st} rx antenna)');

        %% Get Beamforming Weights


        nLayers = 1;
        scOffset = 0;   % no offset
        noRBs = 1;      % average channel conditions over 1 RB to calculate beamforming weights
        [wbs,wue,~] = getBeamformingWeights(hest,nLayers,scOffset,noRBs);
        %% Plot Radiation Patterns
        % Plot the radiation patterns obtained for the UE and the base station.

        % Plot UE radiation pattern
        ueSite.Antenna = clone(channel.ReceiveAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
        ueSite.Antenna.Taper = wue;
        pattern(ueSite,fc,"Size",4);

        % Plot BS radiation pattern
        bsSite.Antenna = clone(channel.TransmitAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
        bsSite.Antenna.Taper = wbs;
        pattern(bsSite,fc,"Size",5);
    end
    fprintf("进度条：%d / %d \n",index,number);
    
   
end

%% 保存数据集

save(filename,'HCSI','userset','bsset','-v7.3');
fprintf("数据保存完毕 \n");
viewer.close();
close all;
clear HCSI userset bsset;

end

%% Local Functions
function [wtx,wrx,D] = getBeamformingWeights(hEst,nLayers,scOffset,noRBs)
    % Get beamforming weights given a channel matrix hEst and the number of
    % layers nLayers. One set of weights is provided for the whole bandwidth.
    % The beamforming weights are calculated using singular value (SVD)
    % decomposition.
    %
    % Only part of the channel estimate is used to get the weights, this is
    % indicated by an offset SCOFFSET (offset from the first subcarrier) and a
    % width in RBs (NORBS).

    % Average channel estimate
    [~,~,R,P] = size(hEst);
    %H = permute(mean(reshape(hEst,[],R,P)),[2 3 1]);

    scNo = scOffset+1;
    hEst = hEst(scNo:scNo+(12*noRBs-1),:,:,:);
    H = permute(mean(reshape(hEst,[],R,P)),[2 3 1]);

    % SVD decomposition
    [U,D,V] = svd(H);
    wtx = V(:,1:nLayers).';
    wrx = U(:,1:nLayers)';
end

function [range] = randpostion(min,max,number)
%RANDPOSTION 此处显示有关此函数的摘要
%   此处显示详细说明

range = (max-min).*rand(number,1) + min;

end


function [H] = dealCSI(hest,RBNum)
    hest_shape = size(hest);
    H = zeros([RBNum,hest_shape(3:end)]);
    hest = mean(hest,2); 
    hest = reshape(hest,[hest_shape(1),hest_shape(3:4)]);
    for m = 1:RBNum
        temp = hest((m-1)*12+1:m*12,:,:);
        H(m,:,:) = mean(temp,1);
    end 
end






