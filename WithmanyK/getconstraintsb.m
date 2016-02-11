function [cl,cu] = getconstraintsb(params)

cl = zeros(params.ncon,1);
cu = zeros(params.ncon,1);	
cl(params.nconeq+1:end) = -0.3+zeros(params.NSU,1);%inf;
cu(params.nconeq+1:end) = 0.3+zeros(params.NSU,1);%inf;           
