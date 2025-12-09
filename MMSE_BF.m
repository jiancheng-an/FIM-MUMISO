function MMSE_BF_output = MMSE_BF(channel, gamma, sigma2)
%% channel: uplink channel N*K
%% gamma: target SINR
%% sigma2: noise power
N = size(channel,1);
K = size(channel,2);
ZF_BF = (channel/(channel'*channel))*sqrt(gamma*sigma2);
lambda = (vecnorm(ZF_BF).^2).';
for mm = 1:500                                                 %%% inner loop
    sigma_matrix = 1/sigma2*diag(lambda);
    d_temp = eye(N)+channel*sigma_matrix*channel';
    for kk = 1:K
        lambda(kk) = real(sigma2/((1+1/gamma)*channel(:,kk)'*inv(d_temp)*channel(:,kk)));
    end
end
sigma_matrix = 1/sigma2*diag(lambda);
d_temp = eye(N)+channel*sigma_matrix*channel';
MMSE_direction = d_temp\channel;
MMSE_direction = MMSE_direction./vecnorm(MMSE_direction);
MM = -abs(channel'*MMSE_direction).^2;
MM(1:(K+1):end) = -MM(1:(K+1):end)/gamma;
power = sigma2*(MM\ones(K,1));
MMSE_BF_output = MMSE_direction*diag(sqrt(power));
end