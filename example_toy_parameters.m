% Demonstrates how to run the synthesis algorithm

% name of the audio waveform
fname = 'speech.wav';

% directory containing the audio waveform
input_directory = pwd;

% directory to save results of the synthesis process
output_directory = [pwd '/example_toy_parameters'];
if ~exist(output_directory,'dir'); 
    mkdir(output_directory);
end

% read parameters
P = synthesis_parameters_toy;

% run synthesis
P.n_iter = 100;
run_spectrotemporal_synthesis(P, fname, input_directory, output_directory);