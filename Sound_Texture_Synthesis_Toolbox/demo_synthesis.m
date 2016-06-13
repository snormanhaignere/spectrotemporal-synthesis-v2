%
% Runs a demo of the McDermott and Simoncelli texture synthesis algorithm.
%
% Uses default parameters and synthesis options set to reduce CPU and
% memory requirements (omits C1 correlations and does not display figures
% at the end).
%
% Parameters are set in synthesis_parameters_fast_demo.m , which also
% specifies the sound file from which the statistics are to be measured;
% you will probably want to change this to your own sound file but this
% demo will work straight away on an example file. 
%

% This code is part of an instantiation of a sound texture synthesis
% algorithm developed with Eero Simoncelli and described in this paper:
%
% McDermott, J.H., Simoncelli, E.P. (2011) Sound texture perception via
% statistics of the auditory periphery: Evidence from sound synthesis.
% Neuron, 71, 926-940. 
%
% Dec 2012 -- Josh McDermott <jhm@mit.edu>

synthesis_parameters_demo;
[synth_sound] = run_synthesis(P);
