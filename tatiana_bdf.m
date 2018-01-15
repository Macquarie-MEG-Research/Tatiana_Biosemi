function tatiana_bdf
%clearvars
%close all
ft_defaults
ft_hastoolbox('fastica', 1);

eegs = dir('*.bdf');
for i =1:length(eegs)
    
    filename = eegs(i).name;
    hdr      = ft_read_header(filename);
    
    %need to replace the labels with something standard
    labels = {'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3',...
        'FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3',...
        'P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz',...
        'Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8',...
        'FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6',...
        'CP4','CP2','P2','P4','P6','P8','P10','P08','PO4','O2'};
    
    cfg                      = [];
    cfg.headerfile           = filename;
    cfg.datafile             = filename;
    cfg.trialfun             = 'ft_trialfun_general';
    cfg.trialdef.triallength = Inf;
    cfg.trialdef.ntrials     = 1;
    cfg                      = ft_definetrial(cfg);
    cfg.continuous           = 'yes';
    cfg.channel              = hdr.label(1:66);
    cfg.reref                = 'yes';
    cfg.refchannel           = {'EXG1','EXG2'}; %reref to linked mastoids
    cfg.baselinewindow       = [-0.2 0]; %FIXME given the task we should make this earlier to avoud CPS messing the baseline
    cfg.lpfilter             = 'yes';
    cfg.lpfreq               = 30;
    cfg.hpfilter             = 'yes';
    cfg.hpfreq               = 1; %FIXME should try something lower to appease reviewers
    cfg.hpfiltord            = 5;
    data                     = ft_preprocessing(cfg);
    
    %Get rid of EXG channels
    cfg         = [];
    cfg.channel = data.label(1:64);
    data        = ft_selectdata(cfg,data);
    
    %Add labels
    data.label(1:64) = labels;
    
    %Downsample by a power of 2 - use 16 to return 128Hz
    cfg            = [];
    cfg.resamplefs = data.fsample/16; % Here we are downsampling from 2048Hz --> 128Hz
    cfg.detrend    = 'yes'; % Helps with low-frequency drift
    data           = ft_resampledata(cfg, data);
    data_clean     = data;
    
    datat = data.trial{1}; %create temp variable with data matrix only
    
    %wICA cleaning - *need wavelet toolbox*
    [wIC,A] = wICA(datat,'fastica',1,0,7); %ica type, mutiplier (threshold),plotting, wavelet cycles, wavelet type
    ahat    = A*wIC;
    nhat    = datat-ahat; %Subtract artifact only from data
    
    %%SEE Castellanos & Makarov, J. Neurosci. Method. 2006 %%
    
    data_clean.trial{1} = nhat;
    clear datat
    
    cfg                     = [];
    cfg.headerfile          = filename;
    cfg.datafile            = filename;
    cfg.trialfun            = 'ft_trialfun_general';
    cfg.trialdef.eventtype  = 'STATUS';
    cfg.trialdef.eventvalue = [211 212 213 215 221 222 223 225];%Trigger numbers
    cfg.trialdef.prestim    = 0.5;
    cfg.trialdef.poststim   = 1;
    cfg                     = ft_definetrial(cfg);
    
    %Here we downsample the matrix representing the trial structure
    cfg.trl(:,1) = round(cfg.trl(:,1)/16); %FIXME this is possibly a bit dirty - is there a better way?
    cfg.trl(:,2) = cfg.trl(:,1)+length(-(cfg.trialdef.prestim):1/data.fsample:cfg.trialdef.poststim)-1;
    cfg.trl(:,3) = round(cfg.trl(:,3)/16);
    %*********************************************
    
    %Epoch the data - first find trial indices corresponding to the events
    %of interest and group
    data       = ft_redefinetrial(cfg,data);
    data_clean = ft_redefinetrial(cfg,data_clean);
    
    cfg = [];
    
    gr_pause_onset_trials = find(ismember(data.trialinfo,211));
    ug_pause_onset_trials = find(ismember(data.trialinfo,221));
    
    gr_pause_start_trials = find(ismember(data.trialinfo,212));
    ug_pause_start_trials = find(ismember(data.trialinfo,222));
    
    gr_pause_end_trials = find(ismember(data.trialinfo,213));
    ug_pause_end_trials = find(ismember(data.trialinfo,223));
    
    gr_pause_null_trials = find(ismember(data.trialinfo,225));
    ug_pause_null_trials = find(ismember(data.trialinfo,215));
    
    fldnm = ['Subject_',num2str(i)]; %Dynamic field name to create struct with all subs in it
    
    cfg.trials                            = gr_pause_onset_trials;
    eeg_raw.(fldnm).gr_pause_onset_data       = ft_redefinetrial(cfg,data);
    eeg_clean.(fldnm).gr_pause_onset_data_clean = ft_redefinetrial(cfg,data_clean);
    
    cfg.trials                            = ug_pause_onset_trials;
    eeg_raw.(fldnm).ug_pause_onset_data       = ft_redefinetrial(cfg,data);
    eeg_clean.(fldnm).ug_pause_onset_data_clean = ft_redefinetrial(cfg,data_clean);
    
    cfg.trials                            = gr_pause_start_trials;
    eeg_raw.(fldnm).gr_pause_start_data       = ft_redefinetrial(cfg,data);
    eeg_clean.(fldnm).gr_pause_start_data_clean = ft_redefinetrial(cfg,data_clean);
    
    cfg.trials                            = ug_pause_start_trials;
    eeg_raw.(fldnm).ug_pause_start_data       = ft_redefinetrial(cfg,data);
    eeg_clean.(fldnm).ug_pause_start_data_clean = ft_redefinetrial(cfg,data_clean);
    
    cfg.trials                          = gr_pause_end_trials;
    eeg_raw.(fldnm).gr_pause_end_data       = ft_redefinetrial(cfg,data);
    eeg_clean.(fldnm).gr_pause_end_data_clean = ft_redefinetrial(cfg,data_clean);
    
    cfg.trials                          = ug_pause_end_trials;
    eeg_raw.(fldnm).ug_pause_end_data       = ft_redefinetrial(cfg,data);
    eeg_clean.(fldnm).ug_pause_end_data_clean = ft_redefinetrial(cfg,data_clean);
    
    cfg.trials                          = gr_pause_null_trials;
    eeg_raw.(fldnm).gr_pause_null_data       = ft_redefinetrial(cfg,data);
    eeg_clean.(fldnm).gr_pause_null_data_clean = ft_redefinetrial(cfg,data_clean);
    
    cfg.trials                          = ug_pause_null_trials;
    eeg_raw.(fldnm).ug_pause_null_data       = ft_redefinetrial(cfg,data);
    eeg_clean.(fldnm).ug_pause_null_data_clean = ft_redefinetrial(cfg,data_clean);
    
    %Do the averaging
    cfg                                 = [];
    eeg_raw_timelock.(fldnm).gr_pause_onset_timelock = ft_timelockanalysis(cfg, eeg_raw.(fldnm).gr_pause_onset_data);
    eeg_raw_timelock.(fldnm).ug_pause_onset_timelock = ft_timelockanalysis(cfg, eeg_raw.(fldnm).ug_pause_onset_data);
    eeg_raw_timelock.(fldnm).gr_pause_start_timelock = ft_timelockanalysis(cfg, eeg_raw.(fldnm).gr_pause_start_data);
    eeg_raw_timelock.(fldnm).ug_pause_start_timelock = ft_timelockanalysis(cfg, eeg_raw.(fldnm).ug_pause_start_data);
    eeg_raw_timelock.(fldnm).gr_pause_end_timelock   = ft_timelockanalysis(cfg, eeg_raw.(fldnm).gr_pause_end_data);
    eeg_raw_timelock.(fldnm).ug_pause_end_timelock   = ft_timelockanalysis(cfg, eeg_raw.(fldnm).ug_pause_end_data);
    eeg_raw_timelock.(fldnm).gr_pause_null_timelock   = ft_timelockanalysis(cfg, eeg_raw.(fldnm).gr_pause_null_data);
    eeg_raw_timelock.(fldnm).ug_pause_null_timelock   = ft_timelockanalysis(cfg, eeg_raw.(fldnm).ug_pause_null_data);
    
    eeg_clean_timelock.(fldnm).gr_pause_onset_timelock_clean = ft_timelockanalysis(cfg, eeg_clean.(fldnm).gr_pause_onset_data_clean);
    eeg_clean_timelock.(fldnm).ug_pause_onset_timelock_clean = ft_timelockanalysis(cfg, eeg_clean.(fldnm).ug_pause_onset_data_clean);
    eeg_clean_timelock.(fldnm).gr_pause_start_timelock_clean = ft_timelockanalysis(cfg, eeg_clean.(fldnm).gr_pause_start_data_clean);
    eeg_clean_timelock.(fldnm).ug_pause_start_timelock_clean = ft_timelockanalysis(cfg, eeg_clean.(fldnm).ug_pause_start_data_clean);
    eeg_clean_timelock.(fldnm).gr_pause_end_timelock_clean   = ft_timelockanalysis(cfg, eeg_clean.(fldnm).gr_pause_end_data_clean);
    eeg_clean_timelock.(fldnm).ug_pause_end_timelock_clean   = ft_timelockanalysis(cfg, eeg_clean.(fldnm).ug_pause_end_data_clean);
    eeg_clean_timelock.(fldnm).gr_pause_null_timelock_clean   = ft_timelockanalysis(cfg, eeg_clean.(fldnm).gr_pause_null_data_clean);
    eeg_clean_timelock.(fldnm).ug_pause_null_timelock_clean   = ft_timelockanalysis(cfg, eeg_clean.(fldnm).ug_pause_null_data_clean);
    
    %**SOME OPTIONAL FIGURES BELOW
    
%     cfg                                      = [];
%     cfg.lpfilter                             = 'yes';
%     cfg.lpfreq                               = 20;
%     %cfg.vartrllength                         = 2;
%     
%     eeg.(fldnm).filt_gr_pause_onset_timelock = ft_preprocessing(cfg,eeg.(fldnm).gr_pause_onset_timelock);
%     eeg.(fldnm).filt_ug_pause_onset_timelock = ft_preprocessing(cfg,eeg.(fldnm).ug_pause_onset_timelock);
%     eeg.(fldnm).filt_gr_pause_start_timelock = ft_preprocessing(cfg, eeg.(fldnm).gr_pause_start_timelock);
%     eeg.(fldnm).filt_ug_pause_start_timelock = ft_preprocessing(cfg, eeg.(fldnm).ug_pause_start_timelock);
%     eeg.(fldnm).filt_gr_pause_end_timelock   = ft_preprocessing(cfg, eeg.(fldnm).gr_pause_end_timelock);
%     eeg.(fldnm).filt_ug_pause_end_timelock   = ft_preprocessing(cfg, eeg.(fldnm).ug_pause_end_timelock);
%     eeg.(fldnm).filt_gr_pause_null_timelock   = ft_preprocessing(cfg, eeg.(fldnm).gr_pause_null_timelock);
%     eeg.(fldnm).filt_ug_pause_null_timelock   = ft_preprocessing(cfg, eeg.(fldnm).ug_pause_null_timelock);
%     
%     eeg.(fldnm).filt_gr_pause_onset_timelock_clean = ft_preprocessing(cfg,eeg.(fldnm).gr_pause_onset_timelock_clean);
%     eeg.(fldnm).filt_ug_pause_onset_timelock_clean = ft_preprocessing(cfg,eeg.(fldnm).ug_pause_onset_timelock_clean);
%     eeg.(fldnm).filt_gr_pause_start_timelock_clean = ft_preprocessing(cfg, eeg.(fldnm).gr_pause_start_timelock_clean);
%     eeg.(fldnm).filt_ug_pause_start_timelock_clean = ft_preprocessing(cfg, eeg.(fldnm).ug_pause_start_timelock_clean);
%     eeg.(fldnm).filt_gr_pause_end_timelock_clean   = ft_preprocessing(cfg, eeg.(fldnm).gr_pause_end_timelock_clean);
%     eeg.(fldnm).filt_ug_pause_end_timelock_clean   = ft_preprocessing(cfg, eeg.(fldnm).ug_pause_end_timelock_clean);
%     eeg.(fldnm).filt_gr_pause_null_timelock_clean   = ft_preprocessing(cfg, eeg.(fldnm).gr_pause_null_timelock_clean);
%     eeg.(fldnm).filt_ug_pause_null_timelock_clean   = ft_preprocessing(cfg, eeg.(fldnm).ug_pause_null_timelock_clean);
    
    % figure
    % cfg        = [];
    % cfg.xlim=[-0.2 1.0];
    % cfg.layout = 'biosemi64.lay';
    % ft_multiplotER(cfg,filt_gr_pause_start_timelock_clean,filt_ug_pause_start_timelock_clean)
    %
    % figure
    % cfg        = [];
    % cfg.xlim=[-0.2 1.0];
    % cfg.layout = 'biosemi64.lay';
    % ft_multiplotER(cfg,filt_gr_pause_onset_timelock,filt_ug_pause_onset_timelock)
    %
    % figure
    % cfg        = [];
    % cfg.xlim=[-0.2 1.0];
    % cfg.layout = 'biosemi64.lay';
    % ft_multiplotER(cfg,filt_gr_pause_onset_timelock,filt_gr_pause_onset_timelock_clean)
    
    % figure
    % cfg        = [];
    % cfg.xlim=[0.12 0.14];
    % cfg.layout = 'biosemi64.lay';
    % ft_topoplotER(cfg,filt_gr_pause_onset_timelock,filt_ug_pause_onset_timelock)
    
%     figure;plot(eeg.(fldnm).filt_ug_pause_onset_timelock.time,mean(eeg.(fldnm).filt_ug_pause_onset_timelock.avg))
%     hold on
%     plot(eeg.(fldnm).filt_ug_pause_onset_timelock_clean.time,mean(eeg.(fldnm).filt_ug_pause_onset_timelock_clean.avg))
%     legend(['Subject ',num2str(i),' dirty'],['Subject ',num2str(i),' clean'])
    
end

%Now we need to get of the pesky .previous field that gets massive. Is this
%a bug or am I not doing the analysis correctly? Perhaps you shouldn't
%overwrite variables?

d=fieldnames(eeg_clean);
for ii=1:length(d)
e=fieldnames(eeg_clean.(d{ii}));
for iii=1:length(e)
f=fieldnames(eeg_clean.(d{ii}).(e{iii}));
eeg_clean.(d{ii}).(e{iii}).cfg.previous=[];
end
end

d=fieldnames(eeg_clean_timelock);
for ii=1:length(d)
e=fieldnames(eeg_clean_timelock.(d{ii}));
for iii=1:length(e)
f=fieldnames(eeg_clean_timelock.(d{ii}).(e{iii}));
eeg_clean_timelock.(d{ii}).(e{iii}).cfg.previous=[];
end
end

d=fieldnames(eeg_raw);
for ii=1:length(d)
e=fieldnames(eeg_raw.(d{ii}));
for iii=1:length(e)
f=fieldnames(eeg_raw.(d{ii}).(e{iii}));
eeg_raw.(d{ii}).(e{iii}).cfg.previous=[];
end
end

d=fieldnames(eeg_raw_timelock);
for ii=1:length(d)
e=fieldnames(eeg_raw_timelock.(d{ii}));
for iii=1:length(e)
f=fieldnames(eeg_raw_timelock.(d{ii}).(e{iii}));
eeg_raw_timelock.(d{ii}).(e{iii}).cfg.previous=[];
end
end

%Save the structs
save eeg_clean eeg_clean
save eeg_raw eeg_raw

save eeg_clean_timelock eeg_clean_timelock
save eeg_raw_timelock eeg_raw_timelock
end


function [wIC,A,W,IC] = wICA(data,varargin)
%--------------- function [wIC,A,W] = wICA(data,varargin) -----------------
%
% Performs ICA on data matrix (row vector) and subsequent wavelet
% thresholding to remove low-amplitude activity from the computed ICs.
% This is useful for extracting artifact-only ICs in EEG (for example), and
% then subtracting the artifact-reconstruction from the original data. 
%
%               >>> INPUTS >>>
% Required: 
%   data = data matrix in row format
% Optional:
%   type = "fastica" or "radical"...two different ICA algorithms based on
%       entropy. "fastica" (default) is parametric, "radical" is nonparametric.
%   mult = threshold multiplier...multiplies the computed threshold from
%       "ddencmp" by this number. Higher thresh multipliers = less
%       "background" (or low amp. signal) is kept in the wICs.
%   plotting = 1 or 0. If 1, plots wIC vs. non-wavelet thresholded ICs
%   Fs = sampling rate, (for plotting...default = 1);
%   L = level set for stationary wavelet transform. Higher levels give
%       better frequency resolution, but less temporal resolution. 
%       Default = 5
%   wavename = wavelet family to use. type "wavenames" to see a list of
%       possible wavelets. (default = "coif5");
%
%               <<< OUTPUTS <<<
%   wIC = wavelet-thresholded ICs
%   A = mixing matrix (inv(W)) (optional)
%   W = demixing matrix (inv(A)) (optional)
%   IC = non-wavelet ICs (optional)
%   
%       * you can reconstruct the artifact-only signals as:
%               artifacts = A*wIC;
%       - upon reconstruction, you can then subtract the artifacts from your
%       original data set to remove artifacts, for instance.
%
% Example:
%  n = rand(10,1000);
%  a = [zeros(1,400),[.5,.8,1,2,2.4,2.5,3.5,5,6.3,6,4,3.2,3,1.7,1,-.6,-2.2,-4,-3.6,-3,-1,0],zeros(1,578)];
%  data = n + linspace(0,2,10)'*a;
%  [wIC,A] = wICA(data,[],5,1);
%  ahat = A*wIC;
%  nhat = data-ahat;
%  err = sum(sqrt((nhat-n).^2));

% By JMS, 11/10/2015
%---------------------------------------------------------------------------------------

% check inputs
if nargin>1 && ~isempty(varargin{1})
type=varargin{1}; else type='fastica';end
if nargin>2 && ~isempty(varargin{2})
mult=varargin{2};else mult=1;end
if nargin>3 && ~isempty(varargin{3})
plotting=varargin{3}; else plotting=0;end
if nargin>4 && ~isempty(varargin{4})
Fs=varargin{4};else Fs=1;end
if nargin>5 && ~isempty(varargin{5})
L=varargin{5}; else L=5;end
if nargin>6 && ~isempty(varargin{6})
wavename=varargin{6}; else wavename='coif5';end

%%FIXM - do some dimensionality reduction?
% run ICA using "fastica" or "radical"
if strcmp(type,'fastica')
    [IC,A,W] = fastica(data,'approach','defl','g','pow3','displayMode','off'); % fastica for parametric...default "pow3" nonlinearity
elseif strcmp(type,'radical')
    [IC,W] = radical(data); % radical ICA for non-parametric
    A = inv(W);
end

% padding data for proper wavelet transform...data must be divisible by
% 2^L, where L = level set for the stationary wavelet transform
modulus = mod(size(data,2),2^L); %2^level (level for wavelet)
if modulus ~=0
    extra = zeros(1,(2^L)-modulus);
else
    extra = [];
end
      
% loop through ICs and perform wavelet thresholding
disp('Performing wavelet thresholding');
for s = 1:size(IC,1)
    if ~isempty(extra)
        sig = [IC(s,:),extra]; % pad with zeros
    else
        sig = IC(s,:);
    end
    [thresh,sorh,~] = ddencmp('den','wv',sig); % get automatic threshold value
    thresh = thresh*mult; % multiply threshold by scalar
    swc = swt(sig,L,wavename); % use stationary wavelet transform (SWT) to wavelet transform the ICs
    Y = wthresh(swc,sorh,thresh); % threshold the wavelet to remove small values
    wIC(s,:) = iswt(Y,wavename); % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
    clear y sig thresh sorh swc 
end

% remove extra padding
if ~isempty(extra)
    wIC = wIC(:,1:end-numel(extra));
end

% plot the ICs vs. wICs
if plotting>0
    disp('Plotting');
    subplot(3,1,1);
        multisignalplot(IC,Fs,'r');
        title('ICs');
    subplot(3,1,2);
        multisignalplot(wIC,Fs,'r');
        title('wICs')
    subplot(3,1,3);
        multisignalplot(IC-wIC,Fs,'r');
        title('Difference (IC - wIC)');
end

end