function tatiana_GM
load('eeg_clean_timelock')
subjects_included = 1:18; % numbers correspond to order in folder

%initialize variable
gr_pause_start_timelock_clean_all = [];
ug_pause_start_timelock_clean_all = [];
gr_pause_onset_timelock_clean_all = [];
ug_pause_onset_timelock_clean_all = [];

%create a string which lists all subjects
for d = 1:length(subjects_included)-1
    gr_pause_start_timelock_clean_all{d} = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d)),'.gr_pause_start_timelock_clean,'];
    ug_pause_start_timelock_clean_all{d} = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d)),'.ug_pause_start_timelock_clean,'];
    gr_pause_onset_timelock_clean_all{d} = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d)),'.gr_pause_onset_timelock_clean,'];
    ug_pause_onset_timelock_clean_all{d} = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d)),'.ug_pause_onset_timelock_clean,'];
    gr_pause_null_timelock_clean_all{d}  = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d)),'.gr_pause_null_timelock_clean,'];
    ug_pause_null_timelock_clean_all{d}  = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d)),'.ug_pause_null_timelock_clean,'];
end

%add one subject after last comma in list
gr_pause_start_timelock_clean_all{d+1} = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d+1)),'.gr_pause_start_timelock_clean'];
ug_pause_start_timelock_clean_all{d+1} = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d+1)),'.ug_pause_start_timelock_clean'];
gr_pause_onset_timelock_clean_all{d+1} = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d+1)),'.gr_pause_onset_timelock_clean'];
ug_pause_onset_timelock_clean_all{d+1} = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d+1)),'.ug_pause_onset_timelock_clean'];
gr_pause_null_timelock_clean_all{d+1}  = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d+1)),'.gr_pause_null_timelock_clean'];
ug_pause_null_timelock_clean_all{d+1}  = ['eeg_clean_timelock.Subject_',num2str(subjects_included(d+1)),'.ug_pause_null_timelock_clean'];

%do the grand average - keep individual to use for statistics
cfg                = [];
cfg.keepindividual = 'yes';%useful for stats?
eval(['gr_pause_start_timelock_clean_all_GM=ft_timelockgrandaverage(cfg,',cell2mat(gr_pause_start_timelock_clean_all),');'])
eval(['ug_pause_start_timelock_clean_all_GM=ft_timelockgrandaverage(cfg,',cell2mat(ug_pause_start_timelock_clean_all),');'])
eval(['gr_pause_onset_timelock_clean_all_GM=ft_timelockgrandaverage(cfg,',cell2mat(gr_pause_onset_timelock_clean_all),');'])
eval(['ug_pause_onset_timelock_clean_all_GM=ft_timelockgrandaverage(cfg,',cell2mat(ug_pause_onset_timelock_clean_all),');'])
eval(['gr_pause_null_timelock_clean_all_GM=ft_timelockgrandaverage(cfg,',cell2mat(gr_pause_null_timelock_clean_all),');'])
eval(['ug_pause_null_timelock_clean_all_GM=ft_timelockgrandaverage(cfg,',cell2mat(ug_pause_null_timelock_clean_all),');'])

%create a .avg field so can do GFP
gr_pause_start_timelock_clean_all_GM.avg = mean(gr_pause_start_timelock_clean_all_GM.individual);
ug_pause_start_timelock_clean_all_GM.avg = mean(ug_pause_start_timelock_clean_all_GM.individual);
gr_pause_onset_timelock_clean_all_GM.avg = mean(gr_pause_onset_timelock_clean_all_GM.individual);
ug_pause_onset_timelock_clean_all_GM.avg = mean(ug_pause_onset_timelock_clean_all_GM.individual);
gr_pause_null_timelock_clean_all_GM.avg  = mean(gr_pause_null_timelock_clean_all_GM.individual);
ug_pause_null_timelock_clean_all_GM.avg  = mean(ug_pause_null_timelock_clean_all_GM.individual);

%calculate GFP - nice for q&d visualisation
cfg                                      = [];
cfg.method                               = 'amplitude';
gr_pause_start_timelock_clean_all_GM_GFP = ft_globalmeanfield(cfg,gr_pause_start_timelock_clean_all_GM);
ug_pause_start_timelock_clean_all_GM_GFP = ft_globalmeanfield(cfg,ug_pause_start_timelock_clean_all_GM);
gr_pause_onset_timelock_clean_all_GM_GFP = ft_globalmeanfield(cfg,gr_pause_onset_timelock_clean_all_GM);
ug_pause_onset_timelock_clean_all_GM_GFP = ft_globalmeanfield(cfg,ug_pause_onset_timelock_clean_all_GM);
gr_pause_null_timelock_clean_all_GM_GFP  = ft_globalmeanfield(cfg,gr_pause_null_timelock_clean_all_GM);
ug_pause_null_timelock_clean_all_GM_GFP  = ft_globalmeanfield(cfg,ug_pause_null_timelock_clean_all_GM);

%plot GFPs
figure;
cfg = [];
ft_singleplotER(cfg,...
    gr_pause_start_timelock_clean_all_GM_GFP,...
    ug_pause_start_timelock_clean_all_GM_GFP,...
    gr_pause_onset_timelock_clean_all_GM_GFP,...
    ug_pause_onset_timelock_clean_all_GM_GFP,...
    gr_pause_null_timelock_clean_all_GM_GFP,...
    ug_pause_null_timelock_clean_all_GM_GFP);

legend({'gr_pause_start_timelock_clean_all_GM_GFP',...
    'ug_pause_start_timelock_clean_all_GM_GFP',...
    'gr_pause_onset_timelock_clean_all_GM_GFP',...
    'ug_pause_onset_timelock_clean_all_GM_GFP',...
    'gr_pause_null_timelock_clean_all_GM_GFP',...
    'ug_pause_null_timelock_clean_all_GM_GFP'},...
    'Interpreter', 'none');

%Plot grand averages with useful comparisons
figure;
cfg        = [];
cfg.xlim   = [-0.5 1.0];
cfg.layout = 'biosemi64.lay';
ft_multiplotER(cfg,gr_pause_start_timelock_clean_all_GM,ug_pause_start_timelock_clean_all_GM)
legend ('Grammatical pause start','Ungrammatical pause start')
figure;
ft_multiplotER(cfg,gr_pause_start_timelock_clean_all_GM,gr_pause_null_timelock_clean_all_GM)
legend ('Grammatical pause start','Grammatical pause null')
figure;
ft_multiplotER(cfg,ug_pause_start_timelock_clean_all_GM,ug_pause_null_timelock_clean_all_GM)
legend ('Ungrammatical pause start','Ungrammatical pause null')
figure;
ft_multiplotER(cfg,gr_pause_onset_timelock_clean_all_GM,ug_pause_onset_timelock_clean_all_GM)
legend ('Grammatical pause onset','Ungrammatical pausen onset')
figure;
ft_multiplotER(cfg,gr_pause_null_timelock_clean_all_GM,ug_pause_null_timelock_clean_all_GM)
legend ('Grammatical pause null','Ungrammatical pause null')
end
