% global lb ub
lb = [-3 0 -3 0]; ub = [0 1.3 0 1.3];
% load('modelA2DnaK_0201'); x = Solution(3,:);
load('pswarmA2Rp_DnaK_06152018')
local = @(v) modelA2DnaK_dose_response_error(x,v);
% [v,fv,~,output]=pso(local,4,[],[],[],[],lb,ub,[],options);
options = optimoptions('particleswarm','SwarmSize',50,'PlotFcn',@pswplotbestf,'Display','iter',...
    'FunctionTolerance',0.01,'MaxStallIterations',8,'MaxIterations',80);

[x, fval] = particleswarm(local, 4,lb,ub,options);

options2 = optimoptions('fmincon','PlotFcn',@optimplotfval,'Display','iter');
Y = fmincon(local,x([15 37 20 22]),[],[],[],[],lb,ub,[],options2)