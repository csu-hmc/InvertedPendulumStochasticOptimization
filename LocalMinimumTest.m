clear all
close all
clc

NSUs = [50 100 150];% [600:100:1000];
obj = zeros(length(NSUs),10);
for i = 1:length(NSUs)
    NSU = NSUs(i);
    for j = 1:10
        filename = ['Solution_', num2str(NSUs(i)), 'SUs', num2str(j), '.mat'];
        load(filename);
        obj(i,j) = result.obj;
    end
end