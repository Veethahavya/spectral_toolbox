function [opt_params] = variogram_opt(y,s)
%VARIOGRAM_OPT finds the optimum variogram parametres for the y.values given an exponential variogram.

%%
npts  = numel(y.x_pts{1});
dx_pts = cell(s.nd,1);
for i=1:s.nd
  dx_pts{i} = y.x_pts{i}*ones(1,npts)-ones(npts,1)*y.x_pts{i}';
end
% h_eff      = sqrt((dx_pts{1}/s.lambda(1)).^2 + (dx_pts{2}/s.lambda(2)).^2);

% Cyy = s.variance * exp(-h_eff);
% Cyy = s.variance * exp(-(sqrt((dx_pts{1}/s.lambda(1)).^2 + (dx_pts{2}/s.lambda(2)).^2)));


%%
init_params = [-20 1 s.lambda(1) s.lambda(2) s.variance];
init_params_lb = [0 0 5 5 100];
init_params_ub = [abs(max(y.values-y.dem))*10 10 150 150 1000];

func = @(init_params)0.5*log(det(init_params(5)*exp(-(sqrt((dx_pts{1}/init_params(3)).^2 + (dx_pts{2}/init_params(4)).^2)))))+0.5*(y.values-(init_params(1)*ones(npts,1) + init_params(2)*y.dem))'*inv(init_params(5)*exp(-(sqrt((dx_pts{1}/init_params(3)).^2 + (dx_pts{2}/init_params(4)).^2))))*(y.values-(init_params(1)*ones(npts,1) + init_params(2)*y.dem));

%options.Algorithm = 'trust-region-reflective';
%options.iterations = 10000;
%options.funcCount = 1000;
%options.FunctionTolerance = 1e-10;
options.Display = 'final';
warning('off')
% end the opt command (next line) with semicolon to supress all output
opt_params = lsqnonlin(func,init_params,[],[],options)

%% Manual iteration to view the obj. fn.
% init_params = [10 -1];
% for i=1:100000
%     mu = init_params(1)*ones(n,1) + init_params(2)*y.dem;
% 
% %     f = (1/sqrt(det(Cyy)))*exp(-0.5*(y.values-mu)'*inv(Cyy)*(y.values-mu));
% %     fn = @(init_params)-log(f);
%     
%     func = @(init_params)0.5*log(det(Cyy))+0.5*(y.values-mu)'*inv(Cyy)*(y.values-mu);
% 
%     fn_val(i) = func(init_params);
%     diff(i) = mean(y.values-mu);
%     
%     if(mean(y.values-mu)<0.0001)
%         opt_params = init_params;
%         break
%     end
%     
%     init_params(1) = init_params(1)-0.1;
%     init_params(2) = init_params(2)+0.01;
% end
% figure(1)
% plot(fn_val)
% figure(2)
% plot(diff)