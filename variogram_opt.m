function [opt_params] = variogram_opt(Cyy,data,DEM,n)
%VARIOGRAM_OPT finds the optimum variogram parametres for the data given an exponential variogram.

%%
init_params = [-20 10];

%mu = init_params(1)*ones(n,1) + init_params(2)*DEM;

% f = (1/sqrt(det(Cyy)))*exp(-0.5*(data-mu)'*inv(Cyy)*(data-mu));
% fn = @(init_params)-log(f);

func = @(init_params)0.5*log(det(Cyy))+0.5*(data-(init_params(1)*ones(n,1) + init_params(2)*DEM))'*inv(Cyy)*(data-(init_params(1)*ones(n,1) + init_params(2)*DEM));

%options = optimoptions(@lsqnonlin,'
% options.Algorithm = 'trust-region-reflective';
% options.iterations = 10000;
% options.funcCount = 1000;
% options.FunctionTolerance = 1e-10;
% options.Display = 'iter-detailed';
opt_params = fminsearch(func,init_params);


%%
% init_params = [10 -1];
% for i=1:100000
%     mu = init_params(1)*ones(n,1) + init_params(2)*DEM;
% 
% %     f = (1/sqrt(det(Cyy)))*exp(-0.5*(data-mu)'*inv(Cyy)*(data-mu));
% %     fn = @(init_params)-log(f);
%     
%     func = @(init_params)0.5*log(det(Cyy))+0.5*(data-mu)'*inv(Cyy)*(data-mu);
% 
%     fn_val(i) = func(init_params);
%     diff(i) = mean(data-mu);
%     
%     if(mean(data-mu)<0.0001)
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