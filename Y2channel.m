function channel = Y2channel(x_n, y_n, z_n, theta, phi, coe_true)
%% x_n: N*1 vector
%% y_n: N*1 vector
%% z_n: N*1 vector
%% theta: L*1 vector
%% phi: L*1 vector
%% coe_true L*K matrix
phase_matrix = x_n*(sin(theta.').*cos(phi.'))+y_n*(sin(theta.').*sin(phi.'))+z_n*cos(theta.');
channel_nor = exp(1i*2*pi*phase_matrix);
channel = channel_nor*coe_true;
end

