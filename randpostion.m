function [range] = randpostion(min,max,number)
%RANDPOSTION 此处显示有关此函数的摘要
%   此处显示详细说明

range = (max-min).*rand(number,1) + min;

end

