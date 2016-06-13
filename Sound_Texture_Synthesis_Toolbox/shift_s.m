%assumes s is a column vector

function new_s = shift_s(s,num_samples)

if num_samples==0
    new_s=s;
elseif num_samples<0
    new_s = [s(-num_samples+1:end); zeros(-num_samples,1)];
    %new_s = [s(-num_samples+1:end); s(1:-num_samples)];
elseif num_samples>0
    new_s = [zeros(num_samples,1); s(1:end-num_samples)];
    %new_s = [s(end-num_samples+1:end); s(1:end-num_samples)];
end