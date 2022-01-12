clc
clear
%%% Matlab code for MIMO-NOMA baesd VLC with NGDPA power allocation
%%% By Chen Chen from EEE NTU

%% Parameter
BW = 10*10^6; % modulation bandwidth: 10 MHz
MI = 0.5; % Modulation index
RP = 0.53; % Responsivity of photo-detector, A/W
I = 2; J =2 ; %16 QAM

R = 2.5;
Ptx = 10; % Average transmitted optical power: 10 W
delta_PD = 0.04; % PD spacing: 4 cm
delta_LED = 1; %Led Spacing: 1 m


for i = 0:0.1:1
    T = round(i*10) + 1;
    norm_offset = i; %Ini bisa kita adjust atau sesuaikan. (0.1 - 1)
    r = norm_offset*R;
    %% System setup
    xt = 0; yt = 0; zt = 3.85;
    xt1 = xt - delta_LED/2; yt1 = yt; zt1 = zt;
    xt2 = xt + delta_LED/2; yt2 = yt; zt2 = zt;

    % Location of user1: center
    xr_u1 = 0; yr_u1 = 0; zr_u1 = 0.85;
    xr1_u1 = xr_u1 - delta_PD/2; yr1_u1 = yr_u1; zr1_u1 = zr_u1; % PD 1
    xr2_u1 = xr_u1 + delta_PD/2; yr2_u1 = yr_u1; zr2_u1 = zr_u1; % PD 2

    % Location of user2
    xr_u2 = r; yr_u2 = 0; zr_u2 = 0.85;
    xr1_u2 = xr_u2 - delta_PD/2; yr1_u2 = yr_u2; zr1_u2 = zr_u2; % PD 1
    xr2_u2 = xr_u2 + delta_PD/2; yr2_u2 = yr_u2; zr2_u2 = zr_u2; % PD 2
    
    % Location of user3
    xr_u3 = xr_u2 + r; yr_u3 = 0; zr_u3 = 0.85;
    xr1_u3 = xr_u3 - delta_PD/2; yr1_u3 = yr_u3; zr1_u3 = zr_u3; % PD 1
    xr2_u3 = xr_u3 + delta_PD/2; yr2_u3 = yr_u3; zr2_u3 = zr_u3; % PD 2

    %% Channel matrix
    h11_u1 = get_channel_DC_gain(xr1_u1,yr1_u1,zr1_u1,xt1,yt1,zt1);
    h12_u1 = get_channel_DC_gain(xr1_u1,yr1_u1,zr1_u1,xt2,yt2,zt2);
    h21_u1 = get_channel_DC_gain(xr2_u1,yr2_u1,zr2_u1,xt1,yt1,zt1);
    h22_u1 = get_channel_DC_gain(xr2_u1,yr2_u1,zr2_u1,xt2,yt2,zt2);

    h11_u2 = get_channel_DC_gain(xr1_u2,yr1_u2,zr1_u2,xt1,yt1,zt1);
    h12_u2 = get_channel_DC_gain(xr1_u2,yr1_u2,zr1_u2,xt2,yt2,zt2);
    h21_u2 = get_channel_DC_gain(xr2_u2,yr2_u2,zr2_u2,xt1,yt1,zt1);
    h22_u2 = get_channel_DC_gain(xr2_u2,yr2_u2,zr2_u2,xt2,yt2,zt2);
    
    h11_u3 = get_channel_DC_gain(xr1_u3,yr1_u3,zr1_u3,xt1,yt1,zt1);
    h12_u3 = get_channel_DC_gain(xr1_u3,yr1_u3,zr1_u3,xt2,yt2,zt2);
    h21_u3 = get_channel_DC_gain(xr2_u3,yr2_u3,zr2_u3,xt1,yt1,zt1);
    h22_u3 = get_channel_DC_gain(xr2_u3,yr2_u3,zr2_u3,xt2,yt2,zt2);
    
    H_u1 = [h11_u1 h12_u1;h21_u1 h22_u1];
    H_u2 = [h11_u2 h12_u2;h21_u2 h22_u2];
    H_u3 = [h11_u3 h12_u3;h21_u3 h22_u3];

    H_u1_inv = inv(H_u1);
    H_u2_inv = inv(H_u2);
    H_u3_inv = inv(H_u3);

    %% Decoding order
    % LED 1
    h1_u1 = h11_u1 + h21_u1;
    h1_u2 = h11_u2 + h21_u2;
    h1_u3 = h11_u3 + h21_u3;
    H1 = [h1_u1 h1_u2 h1_u3];
    [HH1,I1] = sort(H1); % from small to large

    h1_o1 = HH1(3);
    h1_o2 = HH1(2);
    h1_o3 = HH1(1);

    % LED 2
    h2_u1 = h12_u1 + h22_u1;
    h2_u2 = h12_u2 + h22_u2;
    h2_u3 = h12_u3 + h22_u3;
    H2 = [h2_u1 h2_u2 h2_u3];
    [HH2,I2] = sort(H2); % from small to large

    h2_o1 = HH2(3);
    h2_o2 = HH2(2);
    h2_o3 = HH2(1);

    %% Normalized gain difference power allocation (NGDPA)
    % LED 1
    alpha1 = ((h1_o1 - h1_o2 - h1_o3)/h1_o1)^3;
    
    p1_o3 = 1/(1 + alpha1);
    p1_o2 = alpha1*p1_o3;
    p1_o1 = alpha1*p1_o2;
    P1 = [p1_o3 p1_o2 p1_o1];

    p1_u1 = P1(I1(1));
    p1_u2 = P1(I1(2));
    p1_u3 = P1(I1(3));

    % LED 2
    alpha2 = ((h2_o1 - h2_o2 - h2_o3)/h2_o1)^3;

    p2_o3 = 1/(1 + alpha2);
    p2_o2 = alpha2*p2_o3;
    p2_o1 = alpha2*p2_o2;

    if alpha2 == 0
        P2 = [p2_o1 p2_o2 p2_o3];
    else
        P2 = [p2_o3 p2_o2 p2_o1];
    end

    p2_u1 = P2(I2(1));
    p2_u2 = P2(I2(2));
    p2_u3 = P2(I2(3));

    %% AWGN niose
    Prx1_u1 = (H_u1(1,1)+H_u1(1,2))*Ptx;
    Prx2_u1 = (H_u1(2,1)+H_u1(2,2))*Ptx;
    Prx1_u2 = (H_u2(1,1)+H_u2(1,2))*Ptx;
    Prx2_u2 = (H_u2(2,1)+H_u2(2,2))*Ptx;
    Prx1_u3 = (H_u3(1,1)+H_u3(1,2))*Ptx;
    Prx2_u3 = (H_u3(2,1)+H_u3(2,2))*Ptx;

    Pn1_u1 = sigma2_noise(Prx1_u1,BW);
    Pn2_u1 = sigma2_noise(Prx2_u1,BW);
    Pn1_u2 = sigma2_noise(Prx1_u2,BW);
    Pn2_u2 = sigma2_noise(Prx2_u2,BW);
    Pn1_u3 = sigma2_noise(Prx1_u3,BW);
    Pn2_u3 = sigma2_noise(Prx2_u3,BW);

    Pn1_u1_est = ((H_u1_inv(1,1))^2*Pn1_u1 + (H_u1_inv(1,2))^2*Pn2_u1)...
        /(RP*Ptx*MI)^2;

    Pn2_u1_est = ((H_u1_inv(2,1))^2*Pn1_u1 + (H_u1_inv(2,2))^2*Pn2_u1)...
        /(RP*Ptx*MI)^2;

    Pn1_u2_est = ((H_u2_inv(1,1))^2*Pn1_u2 + (H_u2_inv(1,2))^2*Pn2_u2)...
        /(RP*Ptx*MI)^2;

    Pn2_u2_est = ((H_u2_inv(2,1))^2*Pn1_u2 + (H_u2_inv(2,2))^2*Pn2_u2)...
        /(RP*Ptx*MI)^2;

    Pn1_u3_est = ((H_u3_inv(1,1))^2*Pn1_u3 + (H_u3_inv(1,2))^2*Pn2_u3)...
        /(RP*Ptx*MI)^2;

    Pn2_u3_est = ((H_u3_inv(2,1))^2*Pn1_u3 + (H_u3_inv(2,2))^2*Pn2_u3)...
        /(RP*Ptx*MI)^2;

    %% Sum rate (Mb/s)
    if p1_u1 >= p1_u2 && p1_u1 >= p1_u3
        SNR1_u1 = p1_u1/(p1_u2 + p1_u3 + Pn1_u1_est);
        SNR1_u2 = p1_u2/Pn1_u2_est;
        SNR1_u3 = p1_u3/Pn1_u3_est;
    elseif p1_u2 >= p1_u1 && p1_u2 >= p1_u3
        SNR1_u1 = p1_u1/Pn1_u1_est;
        SNR1_u2 = p1_u2/(p1_u1 + p1_u3 + Pn1_u2_est);
        SNR1_u3 = p1_u3/Pn1_u3_est;
    else
        SNR1_u1 = p1_u1/Pn1_u1_est;
        SNR1_u2 = p1_u2/Pn1_u2_est;
        SNR1_u3 = p1_u3/(p1_u1 + p1_u2 + Pn1_u3_est);
    end
    
    if p2_u1 >= p2_u2 && p2_u1 >= p2_u3
        SNR2_u1 = p2_u1/(p2_u2 + p2_u3 + Pn2_u1_est);
        SNR2_u2 = p2_u2/Pn2_u2_est;
        SNR2_u3 = p2_u3/Pn2_u3_est;
    elseif p2_u2 >= p2_u1 && p2_u2 >= p2_u3
        SNR2_u1 = p2_u1/Pn2_u1_est;
        SNR2_u2 = p2_u2/(p2_u1 + p2_u3 + Pn2_u2_est);
        SNR2_u3 = p2_u3/Pn2_u3_est;
    else
        SNR2_u1 = p2_u1/Pn2_u1_est;
        SNR2_u2 = p2_u2/Pn2_u2_est;
        SNR2_u3 = p2_u3/(p2_u1 + p2_u2 + Pn2_u3_est);
    end
    % % Average BER
    % BER1_u1 = theoMQAM(real(SNR1_u1),I,J);
    % BER2_u1 = theoMQAM(real(SNR2_u1),I,J);
    % BER1_u2 = theoMQAM(real(SNR1_u2),I,J);
    % BER2_u2 = theoMQAM(real(SNR2_u2),I,J);
    % 
    % ABER = BER1_u1 + BER2_u1 + BER1_u2 + BER2_u2;

    % Sum rate
    R_NOMA1_u1 = 0.5*BW*log2(1 + SNR1_u1);
    R_NOMA2_u1 = 0.5*BW*log2(1 + SNR2_u1);
    R_NOMA1_u2 = 0.5*BW*log2(1 + SNR1_u2);
    R_NOMA2_u2 = 0.5*BW*log2(1 + SNR2_u2);
    R_NOMA1_u3 = 0.5*BW*log2(1 + SNR1_u3);
    R_NOMA2_u3 = 0.5*BW*log2(1 + SNR2_u3);

    R_NOMA1 = real(R_NOMA1_u1 + R_NOMA1_u2 + R_NOMA1_u3)/10^6; % LED 1
    R_NOMA2 = real(R_NOMA2_u1 + R_NOMA2_u2 + R_NOMA2_u3)/10^6; % LED 2

    R_NOMA(T) = R_NOMA1 + R_NOMA2 % unit: Mbit/s
end

plot(0:0.1:1,R_NOMA)
legend


