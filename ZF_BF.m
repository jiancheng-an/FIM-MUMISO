function ZF_BF_output = ZF_BF(channel, gamma, sigma2)
%% channel: uplink channel N*K
%% gamma: target SINR
%% sigma2: noise power
N = size(channel,1);
K = size(channel,2);
ZF_BF_output = (channel/(channel'*channel))*sqrt(gamma*sigma2);
end