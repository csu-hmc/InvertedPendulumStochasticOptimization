clear all
close all
clc

NSUs = [1 2 3 4 5 6 7 8 9 10 20 30 40 60 70 80 90 110 130 150];
obj = zeros(length(NSUs),10);
for i = 1:length(NSUs)
    NSU = NSUs(i);
    for j = 1:10
        filename = ['Solutions/Solution_', num2str(NSUs(i)), 'SUs', num2str(j), '.mat'];
        load(filename);
        obj(i,j) = result.obj;
    end
end

plot(NSUs,mean(obj,2));