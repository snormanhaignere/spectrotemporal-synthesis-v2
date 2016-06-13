% FUNCTION J_SPECGRAM2(S,SR,SHOW_COLORBAR)
%
% displays spectrogram of signal S with sampling rate SR, showing colorbar
% alongside it if SHOW_COLORBAR is set to 1
%
% uses 20ms windows with 50% overlap, and plots amplitude in dB
%

% Josh McDermott <jhm@mit.edu>

function j_specgram2(s,sr,show_colorbar)

if nargin==2
    show_colorbar=1;
end

win = round(.02*sr);
noverlap = round(win/2);
[S,F,T] = spectrogram(s,win,noverlap,[],sr);
min_t = floor(min(T));
%max_t = floor(max(T));
max_t = max(T);
min_f = floor(min(F));
max_f = floor(max(F));

A = abs(S);
imagesc(T,F,20*log10(A/max(max(A))));
%surf(T,F,20*log10(A/max(max(A))));
%shading flat;view(0,90);
axis xy
invgray = flipud(gray);
colormap(invgray);set(gca,'CLim',[-50 0]);
if show_colorbar
    colorbar
end
set(gca,'XLim', [min_t max_t]);
set(gca,'YLim', [min_f max_f]);