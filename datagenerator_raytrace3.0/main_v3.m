clc; clear; close all;
format long
%% 参数设置
mapfile = "./mapfile/zybmap.osm";
freq = 28e9;  %中心频率
RBNum = 26;   %资源块数
BSTX = [8,4];   %基站天线
UETX = [2,2];    %用户天线
% Bandwidth configuration, required to set the channel sampling rate and for perfect channel estimation
SCS = 30; % subcarrier spacing
NRB = RBNum; % number of resource blocks, 10 MHz bandwidth
%% 基站选点
blocknum = 10; samplenum = 300;
bsPosition_list = zeros(blocknum,2);
bsPosition_list(:,1) = randpostion(32.0516, 32.0293, blocknum); %lat
bsPosition_list(:,2) = randpostion(118.7709, 118.7963, blocknum); %lont
bsHight_list = randpostion(8, 20, blocknum);
bsrange_list = randpostion(200, 200, blocknum);
%% 按block循环
for block_index = 1:blocknum
    filename = "../data_ray3.0/"+"f"+num2str(freq/1e9)+"_NRB"+num2str(RBNum)+"_"+...
    num2str(samplenum)+"_BLOCK"+num2str(block_index)+".mat" %保存路径
    datagenerator_rayt(bsPosition_list(block_index,:),bsHight_list(block_index),...
        bsrange_list(block_index),mapfile,filename,samplenum)
end