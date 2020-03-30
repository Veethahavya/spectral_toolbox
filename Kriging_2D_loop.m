close all
clear all
load heads.mat                                                              % loads head data
load dates.mat                                                              % dates of the run
load estimate_delineator.mat
m = 1;                                                                      % start day
n = 10;                                                                   % end day
figure('Name','Estimate of Study Area');

for i = m:n
    %% definition of prior covariance
    s.model      = 'exponential';                                           % geostatistical model for unknowns
    s.variance   = 528 ;                                                    % geostatistical parameter for exponential model
    s.lambda     = [18.5 18.5];                                             % geostatistical parameter for exponential model
    s.nugget     = 0;                                                       % nugget effect
    %s.kappa      = 1.5;                                                    % shape parameter (for matern covariance only)
    %s.micro      = 0.1;                                                    % microscale smoothing parameter (before nugget)
    s.opt_params = [];
    s.opt = 'true';                                                         % switch to turn variogram optimization on/off
    s.opt_itr = 1;                                                          % number of iterations for variogram optimisation

    %% definition of the grid for the unknowns s
    s.n_pts      = [98  98];                                                % number of unknowns in each direction
    s.d_pts      = [1   1];                                                 % grid spacing in each direction
    s.npts       = prod(s.n_pts);                                           % number of unknowns (not required unless used in problem generation)

    %% definition of mean - _add influence of topo from DEM here?_
    b.model      = 'zero';                                                  % uncertain (this), known, zero or unknown
    b.n          = 1;                                                       % number of base functions
    b.beta_pri   = 0;  %mean for regression fn with dem (array if we have many regression fns)                                                    % prior mean coefficients for trend functions
    b.Qbb        = 1;                                                       % uncertainty of prior mean: covariance matrix for trend coefficients
    b.trend{1}   = ones(s.n_pts);                                           % constant mean function %repeat this line for dem
    b.trends     = b.trend{1}(:);

    %% definition of measurement locations for measurements y
    y.gridtype   = 'irregular';                                             % type of measurement grid (regular or irregular)
    y.npts       = 10;                                                      % number of observations (not required unless used for problem generation)
    col = [94	83	79	98	84	85	97	94	91	84];
    row = [33	71	82	84	43	58	09	17	27	51];
    %obs_names = ['Condors','Wratts','Selmes','Murphys','Pauls','Giffords','OldMCB','MCB','Catch_sh','P_Neal'];
    y.indices = transpose(sub2ind(s.n_pts,row,col));                        % measurement indices in field of unknowns (required for irregular grids) (r-v; c-h)

    %% generation of data set y
    y.error      = 0;                                                       % measurement error (scalar) expressed as variance
    y.values     = transpose(heads(i,:));
%     for j=1:y.npts
%         y.dem(j,1) = (y.values(j,1) + 20 + 10*rand());                  % DEM = GWL + 20 + noise
%     end
    y.dem = [30 7 5 5 25 9 50 50 37 8];
    
    %% kriging method options
    options.superpos = 'fft';                                               % superposition method: fft or standard
    options.solver   = 'standard';                                          % solver method: fft (fft-reg, fft-irreg) or standard
    options.estvar   = 'none';                                              % estimation variance method: full, one-point, speedy or none
    options.plot     = false;                                               % plotting flag: true or false

    %% specific options for fft-based solvers
    options.tol      = 1e-10;                                               % solver relative residual (required for fft-solvers: usually about 1e-10)
    options.maxit    = 200;                                                 % solver maximum iteration number (required for fft-solvers: usually number of measurements or less)
    options.cond     = 1e6;                                                 % solver regularization parameter (required for fft-solvers: usually 1e6)
    options.verbose  = 0;                                                   % solver verbosity (required for fft-solvers: 0, 1 or 2)
    options.kalstr   = 0;                                                   % solver filter strength (required for fft-solvers: zero or some percent of variance)
    options.flag_Strang = 0;                                                % solver method flag (required for fft-solvers: 1 or 0)
    options.maxprime = 7;                                                   % embedding optimization parameter (2,3,5,7,...)
    
    %% calculation of the Estimate
    [estimate,ksi,beta,est_var,s]=general_kriging(s,b,y,options);
    estimate = rot90(estimate,3);
    estimate = flip(estimate,2);
    local_estimate {i} = estimate.*delineator;
    estimateCell{i} = estimate;
    for x=1:y.npts                                                          % to store the corresponding estimates at data locations
        est_at_data(i,x) = estimate(col(x),row(x));
    end
    heads_avg(i) = mean(heads(i,:));                                        % average of heads from data
    
    %% plotting Estimate and the observation points on the estimate plot
    title('Estimated Hydraulic Heads of Study Area')
    pcolor(local_estimate{i})
    set(gca,'Ydir','reverse')
    colorbar
    hold on
    scatter(row,col,10,'red','filled')                                     % plots the observation locations on the pcolor plot
    hold off
    
    %% saving the Local Results into files
    fig1_file_name = sprintf('%s\estimates\estimate_%s.jpg',pwd,dates(i));
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    %savefig(fig_file_name);                                                % saves plot as .fig
    saveas(gcf,fig1_file_name)                                              % saves plot as .jpg
    est_file_name = sprintf('%s\estimates\estimate_%s.mat',pwd,dates(i));
    save(est_file_name,'estimate');
    
end

%% saving Aggregated Results into files
estCell_file_name = sprintf('%s\estimates\estimateCell_%s_to_%s.mat',pwd,dates(m),dates(n));
save(estCell_file_name,'estimateCell');
local_est_file_name = sprintf('%s\estimates\estimates-local_%s_to_%s.mat',pwd,dates(m),dates(n));
save(local_est_file_name,'local_estimate');

%% calculation, plotting, and saving of the Difference in estimates b/w day m and day n - could turn this into Daily diff and make a movie of it
est_diff_switch = 0;
if est_diff_switch == 1
    est_diff = estimateCell{1,m} - estimateCell{1,n};
    est_diff_file_name = sprintf('%s\estimates\estimateDiff_%s_to_%s.mat',pwd,dates(m),dates(n));
    save(est_diff_file_name,'est_diff');
    figure('Name','Estimate Differences');
    title('Estimate Differences')
    ylabel('Head (m)')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    mesh(est_diff)
    fig2_file_name = sprintf('%s\estimates\estimate_difference_%s_to_%s.jpg',pwd,dates(m),dates(n));
    %savefig(fig_file_name);                                                    % saves plot as .fig
    saveas(gcf,fig2_file_name)                                                  % saves plot as .jpg
end

%% calculation, plotting, and saving of Averages over each Location
est_avg_switch = 0;
if est_avg_switch == 1
    for i = m:n
        for j = 1:s.n_pts(1)
            for k = 1:s.n_pts(1)
                est_avg(j,k) = estimateCell{1,i}(j,k);
            end
        end
    end
    est_avg_file_name = sprintf('%s\estimates\estimate_average_%s_to_%s.mat',pwd,dates(m),dates(n));
    save(est_avg_file_name,'est_avg');
    figure('Name','Average of Estimates');
    title('Average of Estimates')
    ylabel('Head (m)')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    mesh(est_avg)
    fig3_file_name = sprintf('%s\estimates\estimate_average_%s_to_%s.jpg',pwd,dates(m),dates(n));
    %savefig(fig_file_name);                                                    % saves plot as .fig
    saveas(gcf,fig3_file_name)                                                  % saves plot as .jpg
end

%% calculation, plotting, and saving of Temporal Averages over the study area
est_avgT_switch = 0;
if est_avgT_switch == 1
    c = 1;
    for i = m:n
        local_est_avg_temporal(c) = mean(local_estimate{1,i}(63:98,:),'all');
        c = c+1;
    end
    est_avgT_file_name = sprintf('%s\estimates\estimateL_average_temporal_%s_to_%s.mat',pwd,dates(m),dates(n));
    save(est_avgT_file_name,'local_est_avg_temporal');
    figure('Name','Temporal Average of Local Estimates');
    title('Temporal Average of Local Estimates')
    ylabel('Head (m)')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    yyaxis left
    plot(m:n,local_est_avg_temporal,'-x')
    yyaxis right
    plot(m:n,heads_avg,'-x')
    set(gca,'xtick',[m:floor((n-m)/15):n],'xticklabel',dates(m:floor((n-m)/15):n))
    fig3_file_name = sprintf('%s\estimates\estimateL_average_temporal_%s_to_%s.jpg',pwd,dates(m),dates(n));
    %savefig(fig_file_name);                                                    % saves plot as .fig
    saveas(gcf,fig3_file_name)                                                  % saves plot as .jpg
end

%% calculation/display of Difference in Data and Estimate
est_data_switch = 0;
if est_data_switch == 1
    fprintf("The Data is:\n")
    disp(heads(m:n,:))
    fprintf("\nThe Estimate at corresponding data locations is:\n")
    disp(est_at_data)
    fprintf("\nThe Absolute Difference in Data and Estimate is:\n")
    disp(abs(heads(m:n,:)-est_at_data))
end