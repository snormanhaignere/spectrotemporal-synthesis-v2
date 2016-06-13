% sets default values of parameters for the McDermott and Simoncelli sound
% texture synthesis algorithm.
%
% The sound file from which the statistics are to be measured is
% specified as a string in the variable P.orig_sound_filename; you will
% probably want to change this to your own sound file.
%
% The folder from which the sound file should be read is specified as a
% string in P.orig_sound_folder.
%
% The folder to which the synthetic result and accompanying files should be
% saved is specified as a string in P.output_folder.
%

% This code is part of an instantiation of a sound texture synthesis
% algorithm developed with Eero Simoncelli and described in this paper:
%
% McDermott, J.H., Simoncelli, E.P. (2011) Sound texture perception via
% statistics of the auditory periphery: Evidence from sound synthesis.
% Neuron, 71, 926-940. 
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>
%
% modified Nov 2013 to allow initialization with arbitrary stored waveform.
% Nov 2013 -- Josh McDermott <jhm@mit.edu>

rand('state',sum(100*clock));
randn('state',sum(100*clock));

%impose stats whose variables are set to 1
P.constraint_set.sub_var = 1; 
P.constraint_set.sub_kurt = 0;
P.constraint_set.env_mean = 1;
P.constraint_set.env_var = 1;
P.constraint_set.env_skew = 1;
P.constraint_set.env_kurt = 0;
P.constraint_set.env_C = 1;
P.constraint_set.env_ac = 0;
P.constraint_set.mod_pow = 1;
P.constraint_set.mod_C1 = 1;
P.constraint_set.mod_C2 = 0;

%impose values of noise stats for those set to 1
P.use_noise_stats.sub_var = 0;
P.use_noise_stats.sub_kurt = 0;
P.use_noise_stats.env_mean = 0;
P.use_noise_stats.env_var = 0;
P.use_noise_stats.env_skew = 0;
P.use_noise_stats.env_kurt = 0;
P.use_noise_stats.env_C = 0;
P.use_noise_stats.env_ac = 0;
P.use_noise_stats.mod_pow = 0;
P.use_noise_stats.mod_C1 = 0;
P.use_noise_stats.mod_C2 = 0;

P.orig_sound_filename = 'Bubbling_water.wav';
P.orig_sound_folder = 'Example_Textures/'; %must be a string, should have a slash at the end. If this is an empty string, Matlab will search its path for the file.
P.output_folder = 'Output_Folder/'; %must be a string, should have a slash at the end. 

P.initialize_with_sound_file = 0; %if 1, initializes synthetic sound to a waveform in the specified file
P.initial_sound_filename = 'white_noise_5s.wav';
P.initial_sound_folder = 'Example_Textures/';

P.write_norm_orig = 1; %save a copy of the original sound file with same rms as synthetic signal
P.write_spec_match_noise = 0; %create spectrally matched noise for comparison 
P.write_repeat=0; %should be 1 to concatenate the result of the synthesis with itself
P.iteration_snapshots=[]; %will save synthetic signal to file at specified iteration numbers

P.display_figures = 1; %controls whether figures are displayed at the end comparing the statistics of original and synthetic sound
P.save_figures = 1; %if 1, figures are saved and closed rather than left on screen

P.max_orig_dur_s = 7;
P.desired_synth_dur_s = 5;
P.measurement_windowing = 2; %1 --> unwindowed (circular boundary handling), 2 --> global window
P.imposition_windowing = 1; %1 --> unwindowed (circular boundary handling), 2 --> global window
P.win_steepness = .5; % must be between 0 and 1; smaller means steeper window edges
P.imposition_method = 1; %which_method=1; %1--> conj grad; 2--> gauss-newton
P.sub_imposition_order = 1; %1--> starts with subband with most power, works out from there; 2--> different starting subband each iteration

P.avg_stat_option = 0; %1 = average stats across examples; 2 = take weighted average of stats for two sounds (for morphs)
P.morph_ratio = .5; % for morphs; this is the weight of the first sound's statistics in the weighted average
P.files_for_stat_avg = {'Writing_with_pen_on_paper.wav','Bubbling_water.wav','Applause_-_enthusiastic2.wav'}; %only first two will be used for morph
P.avg_filename = 'writing_bubbles';

P.num_it=60; %max number of iterations of synthesis loop
P.end_criterion_db=30; %the synthesis process is halted once all stats that are being imposed attain this SNR
P.flag_criterion_db=20; %if a statistic being imposed does not attain this SNR, the log file is flagged
P.omit_converged_stats=0; %1--> statistics are left out of the imposition process once they attain the criterion SNR level
P.omission_criterion_db=35;
P.sig_sub_cutoff_db = 30; %subbands whose power is this amount or more below that of the max subband are left out of SNR calculations
P.incremental_imposition = 0; % for imposition method #2 - this gradually increases a weight on the adjustment
P.incr_imp_it_limit = 20;
P.num_LS_option=1; %method of choosing number of line searches on each iteration of conjugate gradient
P.init_LS_num = 5; %initial number of line searches for conjugate gradient optimization
P.check_imposition_success = 1; %displays before/after error following each iteration of imposition

P.desired_rms = .01; %.1 was too high; some files clipped during wavwrite; .01 prevents clipping but is a little low for laptop speakers
P.audio_sr = 20000;
P.N_audio_channels = 30; %this is the number excluding lowpass and highpass filters on ends of spectrum
P.low_audio_f = 20; %Hz
P.hi_audio_f = 10000; %Hz
P.use_more_audio_filters = 0; % should be 1 if 2x overcomplete, 2 if 4x overcomplete
P.lin_or_log_filters = 1; %1--> log acoustic & mod; 2--> log acoust, lin mod; 3--> lin acoust, log mod; 4--> lin acoust, lin mod

P.env_sr = 400;
P.N_mod_channels = 20; %These next four parameters control the modulation filterbank from which modulation power is measured
P.low_mod_f = 0.5; %Hz
P.hi_mod_f = 200; %Hz
P.use_more_mod_filters=0; % should be 1 if 2x overcomplete,
P.mod_filt_Q_value = 2;
P.use_zp = 0;% 0 means circular convolution; 1 means zeropadding (for modulation filtering)
P.low_mod_f_C12=1; %Hz - this is the lowest frequency in the octave-spaced modulation filterbank used for the C1 and C2 correlations

P.compression_option=1; % should be 0 for no compression, 1 for power compression, 2 for logarithmic compression
P.comp_exponent = .3;
P.log_constant = 10^-12; %this is the constant added prior to taking the log

P.match_env_hist = 0; %set to 1 to match full envelope histogram (in addition to the moments)
P.match_sub_hist = 0; %set to 1 to match full subband histogram (in addition to the moments)
P.n_hist_bins = 128;
P.manual_mean_var_adjustment = 0; %set envelope mean and variance outside of joint optimization procedure

P.C_offsets_to_impose = [1 2 3 5 8 11 16 21]*2^P.use_more_audio_filters;
%P.num_C_offsets = 8;
P.C_offsets_to_impose = [-fliplr(P.C_offsets_to_impose) P.C_offsets_to_impose];

P.mod_C1_offsets_to_impose = [1 2]*2^P.use_more_audio_filters;
P.mod_C1_offsets_to_impose = [-fliplr(P.mod_C1_offsets_to_impose) P.mod_C1_offsets_to_impose];

P.env_ac_intervals_smp = [1 2 3 4 5 6 7 9 11 14 18 22 28 36 45 57 73 92 116 148 187 237 301]; %in samples

%for subband autocorrelation measurement
P.sub_ac_undo_win = 1; % divide by ac of window
P.sub_ac_win_choice = 2; % type of window
P.num_sub_ac_period = 5; %num periods of subband cf over which to match sub ac

P.neg_env_skew = 0; %should be 1 to impose envelope skew of the opposite sign
P.neg_mod_C2 = 0; %should be 1 to flip sign of imaginary part, 2 to flip sign of real part, 3 to flip sign of both


