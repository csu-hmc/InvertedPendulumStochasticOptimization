clear all
close all
clc

NSUs = 50;% [50 100 110 130];% [600:100:1000];
for i = 1:length(NSUs)
    NSU = NSUs(i);
    disp(NSUs(i))
    for j = 1:10
        result = runOPT(NSU);
        filename = ['Solutions/WithConstraint/NoiseCheck_', num2str(NSUs(i)), 'SUs', num2str(j), '.mat'];
        save(filename, 'result')
    end
end
