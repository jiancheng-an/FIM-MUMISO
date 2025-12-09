function [residual_output, label] = SINR_residual(channel, ww, gamma, sigma2)
%% channel: uplink channel N*K
K = size(channel,2);
A = abs(channel'*ww).^2;
for hj = 1:K
    h_sep(hj) = A(hj,hj)*(1/sigma2+1/(gamma*sigma2))-sum(A(hj,:))/sigma2-1;
end
residual_output = sum(h_sep);
label = sum(min(h_sep,0));
end