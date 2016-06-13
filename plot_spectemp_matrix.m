function plot_spectemp_matrix(spectemp_matrix, spec_scales, temp_rates, varargin)

% Plots an image parameterized by spectral/temporal modulation. 
% 
% -- Example --
% temp_mod_rates = logspace(log10(0.5),log10(128),9);
% temp_mod_rates = [temp_mod_rates, -temp_mod_rates];
% spec_mod_scales = logspace(log10(0.125),log10(8),7);
% spectemp_matrix = randn( length(spec_mod_rates), length(temp_mod_rates) );
% plot_spectemp_matrix(spectemp_matrix, spec_mod_scales, temp_mod_rates);

% sort columns by temporal modulation rate
[~,xi] = sort(temp_rates, 'ascend');
temp_rates = temp_rates(xi);
spectemp_matrix = spectemp_matrix(:,xi);

% plot matrix
imagesc(flipud(spectemp_matrix));

% select xticks to mark
if optInputs(varargin, 'XTick')
    
    xticks = varargin{optInputs(varargin, 'XTick') + 1};
    
else
    
    % tick marks for negative values
    xticks = select_int(get(gca, 'XTick'));
    xticks_neg = xticks(temp_rates(xticks) < 0);
    
    % tick marks for corresponding positve values
    xticks_pos = find(ismember(temp_rates, -temp_rates(xticks_neg)));
    
    xticks = [xticks_neg, xticks_pos];
    
end

% mark the chosen positive and negative xticks
set(gca, 'XTick', xticks, 'XTickLabel', temp_rates(xticks));
xlabel('Temp Mod (Hz)');


% integer tick marks
yticks = select_int(get(gca, 'YTick'));
spec_scales_flip = flip(spec_scales);
set(gca, 'YTick', yticks, 'YTickLabel', spec_scales_flip(yticks));
ylabel('Spec Mod (cyc/oct)');


function y = select_int(x)

eps = 1e-6;
is_int = abs(x-round(x)) < eps;
y = x(is_int);
