clc;
clearvars;
close all;
L = 4; %% Number of paths
K = 4; %% Number of users
N = 4; %% Number of elements
sigma2 = 10^(-9.4); %% noise power
gamma = 10^(0.5); %% QoS
C0 = 10^(-61.38/10); %%
N_x = 2;
N_set = 1:N;
x_n = (0.5*mod(N_set-1,N_x)).';
z_n = (0.5*floor((N_set-1)/N_x)).';
rng(42)
theta = rand(L,1)*pi*2/3+pi/6;
phi = rand(L,1)*pi*2/3+pi/6;
angle = rand(K,1)*pi*2;
Distance_X = rand(K,1);
coe = sqrt(1/2)*(randn(L,K) +1i*randn(L,K));
y_n = rand(N,1)*0.5;
location_x = sin(angle).*Distance_X*10;
location_y = cos(angle).*Distance_X*10+20;
location_z = -5;
distance = sqrt(location_x.^2+location_y.^2+location_z^2);
PATH_LOSS = C0*distance.^(-2.2)/L;
coe_true = coe*diag(sqrt(PATH_LOSS));
channel = Y2channel(x_n, y_n, z_n, theta, phi, coe_true);
ww = MMSE_BF(channel, gamma, sigma2);
Range = 1;
gradient = zeros(N,1);
for tt = 1:50
    tic
    %% gradient
    for nn1 = 1:N
        for kk1 = 1:K
            for kk2 = 1:K
                par_1 = x_n(nn1)*(sin(theta.').*cos(phi.'))+z_n(nn1)*cos(theta.');
                par_2 = exp(1i*2*pi*par_1);
                for ll1 = 1:L
                    par_3(ll1) = coe_true(ll1, kk1)'*par_2(ll1)'*ww(nn1, kk2);
                end
                par_4 = 2*pi*(sin(theta.').*sin(phi.'));
                par_5 = channel(:, kk1);
                par_5(nn1, :) = 0;
                par_6 = par_5'*ww(:, kk2);
                for ll1 = 1:L
                    par_7(ll1) = 2*imag(par_6'*par_3(ll1)*exp(-1i*par_4(ll1)*y_n(nn1)))*par_4(ll1);
                    for ll2 = 1:L
                        par_8(ll1, ll2) = imag(par_3(ll1)'*par_3(ll2)*exp(-1i*(par_4(ll2)-par_4(ll1))*y_n(nn1)))*(par_4(ll2)-par_4(ll1));
                    end
                end
                gradient_temp(kk1,kk2) = sum(par_7)+sum(sum(par_8));
            end
        end
        gradient(nn1) = sum(diag(gradient_temp))*(1+gamma)-gamma*sum(sum(gradient_temp));
    end
    gradient = gradient./norm(gradient)*Range*N;
    %%
    [residual(1), ~] = SINR_residual(channel, ww, gamma, sigma2);
    power_temp(1) = norm(ww)^2;
    step_set = 0:0.001:1;
    for iii = 2:length(step_set)
        y_n_temp = y_n + gradient*step_set(iii);
        y_n_temp = min(max(y_n_temp,0),Range);
        channel = Y2channel(x_n, y_n_temp, z_n, theta, phi, coe_true);
        [residual(iii), label] = SINR_residual(channel,ww, gamma, sigma2);
        ww1 = MMSE_BF(channel, gamma, sigma2);
        power_temp(iii) = norm(ww1)^2;
        if label < -30 || residual(iii) <= residual(iii-1) || power_temp(iii) >= power_temp(iii-1)
            break
        end
    end
    [~,d_z] = min(power_temp);
    step = step_set(d_z);
    y_n = step * gradient + y_n;
    y_n = min(max(y_n,0),Range);
    channel = Y2channel(x_n, y_n, z_n, theta, phi, coe_true);
    %% MU beamforming
    ww = MMSE_BF(channel, gamma, sigma2);
    power_T(tt) = norm(ww)^2;
    figure(1)
    p1 = plot(tt, y_n(1),'ko','markersize',5,'linewidth',1.5,'MarkerFaceColor','c');
    hold on
    p2 = plot(tt, y_n(2),'rh','markersize',4.25,'linewidth',1.5,'MarkerFaceColor','c');
    hold on
    p3 = plot(tt, y_n(3),'bp','markersize',3.5,'linewidth',1.5,'MarkerFaceColor','c');
    hold on
    par_4 = plot(tt, y_n(4),'ms','markersize',5.5,'linewidth',1.5,'MarkerFaceColor','c');
    hold on
    toc
end
figure(1)
plot([0 50],[1 1],'g--','linewidth',1.5)
legend([p1 p2 p3 par_4], 'Element 1', 'Element 2', 'Element 3', 'Element 4','location','northeast');
x = [0.25,0.2];
y = [0.85,0.9];
a = annotation('textarrow',x,y,'String',' y_m_a_x = \lambda','fontsize',16,'linewidth',1.5,'linestyle','-.');
xlabel('Number of iterations')
ylabel('Position [\times\lambda]')
axis([0 50 0 1])
set(gca,'fontsize',16)

figure(2);
plot(10*log10(power_T))