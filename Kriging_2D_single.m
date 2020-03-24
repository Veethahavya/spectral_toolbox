clear all
load heads.mat                                                          % loads head data
load dates.mat                                                          % dates of the run

%% definition of prior covariance
s.model      = 'exponential';                                           % geostatistical model for unknowns
s.variance   = 528;                                                     % geostatistical parameter for exponential model
s.lambda     = [18.5 18.5];                                             % geostatistical parameter for exponential model
s.nugget     = 0;                                                       % nugget effect
%s.kappa     = 1.5;                                                     % shape parameter (for matern covariance only)
%s.micro     = 0.1;                                                     % microscale smoothing parameter (before nugget)
s.opt_params = [];
s.opt = 'false';                                                         % switch to turn variogram optimization on/off
s.opt_itr = 100;                                                        % number of iterations for variogram optimisation

%% definition of the grid for the unknowns s
s.n_pts      = [98  98];                                                % number of unknowns in each direction
s.d_pts      = [1   1];                                                 % grid spacing in each direction
s.npts       = prod(s.n_pts);                                           % number of unknowns (not required unless used in problem generation)

%% definition of mean
b.model      = 'zero';                                                  % uncertain, known, zero or unknown
b.n          = 1;                                                       % number of base functions
b.beta_pri   = 0;                                                       % prior mean coefficients for trend functions
b.Qbb        = 1;                                                       % uncertainty of prior mean: covariance matrix for trend coefficients
b.trend{1}   = ones(s.n_pts);                                           % constant mean function
b.trends     = b.trend{1}(:);

%% definition of measurement locations for measurements y
y.gridtype   = 'irregular';                                             % type of measurement grid (regular or irregular)
y.npts       = 10;                                                      % number of observations (not required unless used for problem generation)
col = [94 83 79 98 84 85 97 94 91 84];
row = [33 71 82 84 43 58 09 17 27 51];
y.indices    = transpose(sub2ind(s.n_pts,row,col));                     % measurement indices in field of unknowns (required for irregular grids) (r-v; c-h)

%% definition of data set y
y.error      = 0;                                                       % measurement error (scalar) expressed as variance
y.values     = transpose(heads(1,:));
% y.edk_dem = y.values + 20;
for i=1:y.npts
    y.edk_dem(i,1) = (y.values(i,1) + 10*rand());
end

%% kriging method options
options.superpos = 'fft';                                               % superposition method: fft or standard
options.solver   = 'standard';                                          % solver method: fft (fft-reg, fft-irreg) or standard
options.estvar   = 'none';                                              % estimation variance method: full, one-point, speedy or none
options.plot     = true;                                                % plotting flag: true or false

%% specific options for fft-based solvers
options.tol      = 1e-10;                                               % solver relative residual (required for fft-solvers: usually about 1e-10)
options.maxit    = 200;                                                 % solver maximum iteration number (required for fft-solvers: usually number of measurements or less)
options.cond     = 1e6;                                                 % solver regularization parameter (required for fft-solvers: usually 1e6)
options.verbose  = 0;                                                   % solver verbosity (required for fft-solvers: 0, 1 or 2)
options.kalstr   = 0;                                                   % solver filter strength (required for fft-solvers: zero or some percent of variance)
options.flag_Strang = 0;                                                % solver method flag (required for fft-solvers: 1 or 0)
options.maxprime = 7;                                                   % embedding optimization parameter (2,3,5,7,...)

%% calculation of the Estimate
[estimate,ksi,beta,est_var,s] = general_kriging(s,b,y,options);
estimate = rot90(estimate,3);
estimate = flip(estimate,2);
load estimate_delineator.mat
local_estimate = estimate.*delineator;

%% plotting the observation points on the estimate plot
% [row,col] = ind2sub(s.n_pts,y.indices);
title('Estimated Hydraulic Heads of Study Area')
pcolor(local_estimate)
set(gca,'Ydir','reverse')
colorbar
hold on
scatter(row,col,10,'red','filled')                                     % plots the observation locations on the pcolor plot
hold off

%% calculation of Difference in Data and Estimate
est_data_switch = 0;
if est_data_switch == 1
    fprintf("The Data is:\n")
    disp(y.values')
    fprintf("\nThe Estimate at corresponding data locations is:\n")
    for x=1:y.npts
        est_at_data(x) = estimate(col(x),row(x));
    end
    disp(est_at_data)
    fprintf("\nThe Absolute Difference in Data and Estimate is:\n")
    disp(abs(y.values'-est_at_data))
    figure(2)
    plot((y.values'-est_at_data),'-x')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end

%% LOO Error calculation
loo_switch = 0;
if loo_switch == 1
    for i = 1:10
       y.npts = 9;
       loo_row = row(1);
       loo_col = col(1);
       loo_val = y.values(1);

       row_tmp = row(2:10);
       col_tmp = col(2:10);
       y.indices = transpose(sub2ind(s.n_pts,row_tmp,col_tmp));
       y.values = y.values(2:10);

       loo_estimate = general_kriging(s,b,y,options);
       loo_estimate = rot90(loo_estimate,3);
       loo_estimate = flip(loo_estimate,2);
       loo_error(i) = loo_estimate(loo_col,loo_row) - loo_val;

       row = circshift(row,1);
       col = circshift(col,1);
       y.values = transpose(heads(1,:));
       y.values = circshift(y.values,i);
 % The observation "Old_MCB" at (97,09) has very high impact on the estimate
    end
    plot(loo_error)
end

%% Comparision of estimate with initial_heads from modflow
init_h_comp_switch = 0;
    if init_h_comp_switch == 1
        load initial_heads
        surf(local_estimate-initial_heads)
        colorbar
        set(gca,'Ydir','reverse')
        hold on
        scatter(row,col,10,'red','filled')
        hold off
        mean_diff_init_and_est = mean(abs(local_estimate(63:98,:)-initial_heads(63:98,:)),'all');
    end
    
%% plotting of the optimal params through iterations
plt_Sopt_switch = 1;
if plt_Sopt_switch == 1 && isequal(s.opt,'true')
    plt_betas = cell2mat(s.opt_params);
    plot(plt_betas(1:2:200))
    hold on
    plot(plt_betas(2:2:200))
    legend('beta1','beta2')
end