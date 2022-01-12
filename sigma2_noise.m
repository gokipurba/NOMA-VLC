
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function for the calculation of total noise power %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function noise_power = sigma2_noise(Prx,symbol_rate)

% Parameters
q = 1.6e-19; % Electronic c harge
R = 0.53; % Responsivity of PD, A/W

% I_bg = 5100e-6; % Background current: direct sun light without optical filter, A
% I_bg = 1000e-6; % Background current: direct sun light with optical filter, A
% I_bg = 190e-6; % Background current: indirect sun light with optical filter, A
I_bg = 58e-6; % Background current: incandescent light and fluorescent light with optical filter, A

I2 = 0.562; % Noise bandwidth factor
B = symbol_rate; % Equivalent noise bandwidth
k = 1.38e-23; % Boltzmann constant
T = 300; % Kelvin, K
G = 10; % Open loop voltage gain
eta = 112e-12 / 1e-4; % Fixed capacitance, F/m
area_PD = 10^(-4); % Physical area of the PD, 1 cm^2
Gamma = 1.5; % Field-effect transistor (FET) channel noise factor
g_m = 30e-3; % FET transconductance, Siemens
I3 = 0.0868;

%% Total noise power calculation
% Shot noise power
sigma2_shot = 2*q*(R*Prx + I_bg*I2)* B;

% Thermal noise power
sigma2_thermal = 8*pi*k*T*eta*area_PD*B^2*(I2/G + 2*pi*Gamma*eta*area_PD*I3*B/g_m);

% Total noise power
noise_power = sigma2_shot + sigma2_thermal;

