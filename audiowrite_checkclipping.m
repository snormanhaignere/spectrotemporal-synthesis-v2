function audiowrite_checkclipping(fname, wav, sr, varargin)

% audiowrite_checkclipping(fname, wav, sr, varargin)
% 
% Wrapper for audiowrite that checks the waveform will not clip
% 
% -- Example --
% 
% fname = ['temp'  num2str(randi(1e9), '%9d') '.wav'];
% sr = 20e3;
% wav = randn(sr,1);
% audiowrite_checkclipping(fname, wav/10, sr);
% delete(fname);
% audiowrite_checkclipping(fname, wav, sr);

if any(abs(wav(:))>1)
    error('Clipping in attempt to write file: \n%s\n\n', fname);
end
audiowrite(fname, wav, sr, varargin{:});