function [] = figure_test()
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
figure
fprintf('123\n\n');
fplot(@(x) x.^2,[-10 10])
end

