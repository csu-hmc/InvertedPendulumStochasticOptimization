clear all
close all
clc

optimizer = 'IPOPT';
NSUs = [9:1:10 20:10:100]; %[50:10:100];%
for i = 1:length(NSUs)
    NSU = NSUs(i);
    disp(NSUs(i))
    for j = 1:10
        [result] = runOpt(NSU, optimizer);
        if strcmp(result(end).params.optimizer, 'IPOPT')
            while result(end).info ~= 0
                rng('shuffle')
                [result] = runOpt(NSU, optimizer);
            end
        elseif strcmp(result(end).params.optimizer, 'SNOPT')
            while result(end).info ~= 1
                rng('shuffle')
                [result] = runOpt(NSU, optimizer);
            end
        else
            error('No other optimizer')
        end
        filename = ['Solutions/Solution_', num2str(NSUs(i)), 'SUs', num2str(j), '.mat'];
        save(filename, 'result')
    end
end