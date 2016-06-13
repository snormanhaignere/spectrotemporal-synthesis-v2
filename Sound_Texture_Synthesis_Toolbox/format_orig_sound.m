function [orig_sound] = format_orig_sound(P, file_n)
% [ORIG_SOUND] = FORMAT_ORIG_SOUND(P, FILE_N)
%
% formats original sound file for texture statistic measurement
%
% parameters P can be set with SYNTHESIS_PARAMETERS.m
%
% in addition, P.ORIG_SOUND_FILENAME or P.FILES_FOR_STAT_AVG must be
% specified (it is not set by SYNTHESIS_PARAMETERS.m)
%
% FILE_N is used when averaging the statistic of multiple sounds, to denote
% which sound out of the list in P.FILES_FOR_STAT_AVG needs to be
% formatted.
%

% This code is part of an instantiation of a sound texture synthesis
% algorithm developed with Eero Simoncelli and described in this paper:
%
% McDermott, J.H., Simoncelli, E.P. (2011) Sound texture perception via
% statistics of the auditory periphery: Evidence from sound synthesis.
% Neuron, 71, 926-940. 
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>

if P.avg_stat_option==0
    [temp,sr] = wavread([P.orig_sound_folder P.orig_sound_filename]);
elseif P.avg_stat_option>=1
    [temp,sr] = wavread([P.orig_sound_folder P.files_for_stat_avg{file_n}]);
end
if size(temp,2)==2
    temp = temp(:,1); %turn stereo files into mono
end
if sr ~= P.audio_sr
    temp = resample(temp, P.audio_sr, sr);
end
if rem(length(temp),2)==1
    temp = [temp; 0];
end

ds_factor=P.audio_sr/P.env_sr;

if P.max_orig_dur_s==-1 ||  P.max_orig_dur_s > length(temp)/P.audio_sr %use full signal to measure stats
    new_l = floor(length(temp)/ds_factor/2)*ds_factor*2; %to accomodate downsampling of envelope
else
    new_l = ceil(P.max_orig_dur_s*P.audio_sr/ds_factor/2)*ds_factor*2; %to accomodate downsampling of envelope
    if new_l > length(temp)
        new_l = floor(P.max_orig_dur_s*P.audio_sr/ds_factor/2)*ds_factor*2;
    end
end
orig_sound = temp(1:new_l);

orig_sound = orig_sound/rms(orig_sound)*P.desired_rms;
