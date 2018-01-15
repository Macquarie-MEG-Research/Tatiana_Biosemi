function tatiana_stats

layout='/Users/mq20096022/Downloads/fieldtrip-master/template/layout/biosemi64.lay';
neighbours=load('/Users/mq20096022/Documents/Matlab/fieldtrip-master/template/neighbours/biosemi64_neighb.mat');

cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to
% evaluate the effect at the sample level
cfg.latency = [0.08 .12]; %time window of interest - here its the N1
cfg.avgovertime = 'yes';
cfg.parameter = 'individual'; %Because we're using the GM struct as input
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.01;         % alpha level of the sample-specific test statistic that
% will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the
% permutation distribution.
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
% required for a selected sample to be included
% in the clustering algorithm (default=0).
cfg.neighbours = neighbours.neighbours;   % see below
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.correcttail = 'prob';% alpha level of the permutation test corrects by multiplying the p-valuses by 2. Makes interpretation easier
cfg.numrandomization = 500;      % number of draws from the permutation distribution

%Design matrix
design = zeros(1,size(gr_pause_null_timelock_clean_all_GM.individual,1) + size(gr_pause_null_timelock_clean_all_GM.individual,1));
design(1,1:size(gr_pause_null_timelock_clean_all_GM.individual,1)) = 1;
design(1,(size(gr_pause_null_timelock_clean_all_GM.individual,1)+1):(size(gr_pause_null_timelock_clean_all_GM.individual,1) + size(gr_pause_null_timelock_clean_all_GM.individual,1)))= 2;

cfg.design = design;             % design matrix
cfg.ivar  = 1;                   % number or list with indices indicating the independent variable(s)

[stat] = ft_timelockstatistics(cfg, gr_pause_onset_timelock_clean_all_GM,ug_pause_onset_timelock_clean_all_GM);

ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path - better colour map than jet

cfg = [];
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.layout = layout;
cfg.subplotsize = [1 1]; % because we're averaging over time we can set this to 1
ft_clusterplot(cfg, stat);
colormap(flipud(brewermap(64, 'RdBu')))

end
