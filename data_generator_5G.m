clc;clear;
%% parameter setting
number = input("Please input the number of samples: \n"); %样本数量
freq = 2.1e9;  
RBNum = 52;   
BSTX = [8,8];   %BS Anttena
UETX = [2,2];    %UE Anttena
%% 基站和用户坐标区域设置
format long
bsPosition = [32.118980, 118.928961]; % Lat, lon
bsHight =10;                       %bs hight(m)
bsAntSize = BSTX;                   
bsArrayOrientation = [0 0].';       
lat1 = 32.119000; lat2 = 32.119950;     
lont1 = 118.928017; lont2 = 118.930270; 
mapfile = "nyxl.osm"       %choose the map file
filename = "./data_ray/"+"f"+num2str(freq/1e9)+"_NRB"+num2str(RBNum)+"_"+num2str(number)+"_BLOCK1.mat" %path to save dataset

%% 随机
HCSI = zeros(number,RBNum,64,4);
lat = randpostion(lat1,lat2,number);
lont = randpostion(lont1,lont2,number);
UEhight = randpostion(0.5,1.5,number);
UEazi = randpostion(0,360,number);
UEele = randpostion(0,90,number);

%% 循环生成CSI
for index = 1:number

    fc = freq;                           
    uePosition = [lat(index),lont(index)];  
    ueAntSize = UETX;                   
    ueArrayOrientation = [UEazi(index) UEele(index)].';      
    reflectionsOrder = 1;                 

    SCS = 15; 
    NRB = RBNum;

    if exist('viewer','var') && isvalid(viewer) 
        viewer.clearMap();
    else
        viewer = siteviewer("Basemap","openstreetmap","Buildings",mapfile); 
    end



    bsSite = txsite("Name","Base station", ...
        "Latitude",bsPosition(1),"Longitude",bsPosition(2),...
        "AntennaAngle",bsArrayOrientation(1:2),...
        "AntennaHeight",bsHight,...  
        "TransmitterFrequency",fc);

    ueSite = rxsite("Name","UE", ...
        "Latitude",uePosition(1),"Longitude",uePosition(2),...
        "AntennaHeight",UEhight(index),... 
        "AntennaAngle",ueArrayOrientation(1:2));

    
    bsSite.show();
    ueSite.show();
   

    rays = raytrace(bsSite,ueSite,"NumReflections",0:reflectionsOrder,"Type","pathloss");
    
    if isempty(rays{1})
        continue
    else
        plot(rays{1})


        pathToAs = [rays{1}.PropagationDelay]-min([rays{1}.PropagationDelay]);  
        avgPathGains  = -[rays{1}.PathLoss];                                 
        pathAoDs = [rays{1}.AngleOfDeparture];       % AoD of each ray
        pathAoAs = [rays{1}.AngleOfArrival];         % AoA of each ray
        isLOS = any([rays{1}.LineOfSight]);                                     


        channel = nrCDLChannel;
        channel.DelayProfile = 'Custom';
        channel.PathDelays = pathToAs;
        channel.AveragePathGains = avgPathGains;
        channel.AnglesAoD = pathAoDs(1,:);      
        channel.AnglesZoD = 90-pathAoDs(2,:);   
        channel.AnglesAoA = pathAoAs(1,:);       
        channel.AnglesZoA = 90-pathAoAs(2,:);    
        channel.HasLOSCluster = isLOS;
        channel.CarrierFrequency = fc;
        channel.NormalizeChannelOutputs = false; 
        channel.NormalizePathGains = false;    
        c = physconst('LightSpeed');
        lambda = c/fc;
        ueArray = phased.URA('Size',ueAntSize(1:2),'ElementSpacing', 0.5*lambda*[1 1]);            
        channel.ReceiveAntennaArray = ueArray;
        channel.ReceiveArrayOrientation = [ueArrayOrientation(1); (-1)*ueArrayOrientation(2); 0]; 
        bsArray = phased.URA('Size',bsAntSize(1:2),'ElementSpacing', 0.5*lambda*[1 1],'Element',phased.NRAntennaElement);
        channel.TransmitAntennaArray = bsArray;
        channel.TransmitArrayOrientation = [bsArrayOrientation(1); (-1)*bsArrayOrientation(2); 0];   
        ofdmInfo = nrOFDMInfo(NRB,SCS);
        channel.SampleRate = ofdmInfo.SampleRate;

        channel.ChannelFiltering = false;
        [pathGains,sampleTimes] = channel();
        pg=permute(pathGains,[2 1 3 4]);
        if isLOS
            pg = [sum(pg(1:2,:,:,:)); pg(3:end,:,:,:)];
        end
        pg = abs(pg).^2;
        pathFilters = getPathFilters(channel);
        nSlot = 0;
        [offset,~] = nrPerfectTimingEstimate(pathGains,pathFilters);
        hest = nrPerfectChannelEstimate(pathGains,pathFilters,NRB,SCS,nSlot,offset,sampleTimes);

        HCSI(index,:,:,:) = dealCSI(hest,RBNum,BSTX,UETX);

        nLayers = 1;
        scOffset = 0; 
        noRBs = 1;      
        [wbs,wue,~] = getBeamformingWeights(hest,nLayers,scOffset,noRBs);
      
        ueSite.Antenna = clone(channel.ReceiveAntennaArray);
        ueSite.Antenna.Taper = wue;
        pattern(ueSite,fc,"Size",4);

        bsSite.Antenna = clone(channel.TransmitAntennaArray);
        bsSite.Antenna.Taper = wbs;
        pattern(bsSite,fc,"Size",5);
    end
    X = sprintf("进度条：%d / %d",index,number);
    disp(X)
   
end

%% save the dataset

save(filename,'HCSI','lat','lont','-v7.3')
X2 = sprintf("The dataset is saved successfully!");
disp(X2)

%% Local Functions
function [wtx,wrx,D] = getBeamformingWeights(hEst,nLayers,scOffset,noRBs)

    [~,~,R,P] = size(hEst);
    scNo = scOffset+1;
    hEst = hEst(scNo:scNo+(12*noRBs-1),:,:,:);
    H = permute(mean(reshape(hEst,[],R,P)),[2 3 1]);
    [U,D,V] = svd(H);
    wtx = V(:,1:nLayers).';
    wrx = U(:,1:nLayers)';
end

function [range] = randpostion(min,max,number)

range = (max-min).*rand(number,1) + min;
end

function [H] = dealCSI(hest,RBNum,BSTX,UETX)
    hest = permute(hest,[1,2,4,3]); 
    H = zeros(RBNum,prod(BSTX),prod(UETX));
    hest = mean(hest,2); 
    hest = reshape(hest,[RBNum*12,prod(BSTX),prod(UETX)]);
    for m = 1:RBNum
        temp = hest((m-1)*12+1:m*12,:,:);
        H(m,:,:) = mean(temp,1);
    end 
end



