clear all
close all
clc

NSUs = [1:1:10 20:10:100]; %
for i = 1:length(NSUs)
    NSU = NSUs(i);
    disp(NSUs(i))
    for j = 1:10
        [~,result] = runOPT(NSU);
        while result(7).info ~= 0
            rng('shuffle')
            [~,result] = runOpt(NSU);
        end
        filename = ['Solutions/LinFeed/Solution_', num2str(NSUs(i)), 'SUs', num2str(j), '.mat'];
        save(filename, 'result')
    end
end
