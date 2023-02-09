%% Mass of our Aircraft

taper_ratio = 0.55; % Lambda
S_ht = 2.31 * 10.7639; % Horizontal tail area in ft^2
S_vt = 1.48 * 10.7639; % Vertical tail area in ft^2
V_stall = 55; % Stall speed in knots (from FAA Performance google sheet)
CLalpha = 5.4113; % lift-curve slope for NACA 644-421 (1/rad)
S = 11.6 * 10.7639; % Wing area in ft^2
AR = 10; % Aspect Ratio of wing
rho = 0.0023769; % sea level density in slug/ft^3
g = 32.17; % Earth gravity in ft/s^2
b = (AR*S)^(1/2); % span of wing

Vx = 70; % knots (from FAA Performance google sheet)
Vy = 80; % knots (from FAA Performance google sheet)
V_cruise = Vy + (Vy - Vx);

Q = 0.5*rho*(V_cruise)^2; % dynamic pressure
t_c = 0.21; % percentage of wing thickness to chord ratio
Wo = 2241; % design gross weight (lbf)
SF = 1.5; % Safety factor
LF = 3.8*g; % Load factor (SR-22) (in g's)
nz = SF * LF; % ultimate load factor
AR_vt = 1.84; % Vertical tail aspect ratio
taper_vt = 0.55; % Vertical tail taper ratio
Lambda = 4.5; % wing sweep angle at 25% MGC (deg)
Lambda_ht = 4.5; % horizontal tail sweep at 25% MGC (deg) (same as wing)
Lambda_vt = 18; % vertical tail sweep at 25% MGC (deg) (same as wing)
nl = nz; % ultimate landing load factor (SR-22) *ASSUMED = nz*
Wl = 0.7*Wo; % design landing weight in lbf
l_FS = 10.1667; % length of fuselage structure (forward bulkhead to aft frame) in ft (SR-22)
l_HT = 12.4145; % horizontal tail arm, from wing C/4 to HT C/4 in ft
d_FS = 4.167; % depth/height of fuselage structure in ft (SR-22) 
deltaP = 8; % cabin pressure differential, in psi (typically 8psi)
V_P = 137; % Volume of pressurized cabin section in ft^3 (SR-22)
N_occ= 2; % number of occupants
M = 0.21; % mach number at cruise
D = 48.75; % Max Fuselage diameter (in)
L1 = 5.50; %length of nose section in ft(SR22)
L2 =8.22; %Length of the center section in ft(SR22)
L3 =8.70; %Length of nose section in ft (SR22)
N_Eng = 1; % number of engines on airplane
W_Eng = 800 * 2.205; % weight of each installed engine in lbf
Ftail = 0; % =0 for conventional tail, =1 for T-tail
M_h= 32.6192; % mass of hydrogen
M_o= 11.02; % mass of oxygen
M_ot= 22.04; % mass of oxygen tank
M_ht= 52.896; % mass of hydrogen tank

% Wing Weight
Wing_WT = 0.036*S^(0.758)*(AR/(cosd(Lambda)^2))^(0.6)*Q^(0.006)*taper_ratio^(0.04)*(100*t_c/(cosd(Lambda)))^(-0.3)*(nz*Wo)^(0.49);


% Horizontal Tail (HT) Weight
W_HT = 0.016*(nz*Wo)^0.414*Q^0.168*S_ht^(0.896)*((100*t_c)/cosd(Lambda_ht))^(-0.12);


% Vertical Tail (VT) Weight
W_VT = 0.073 * (1+0.2*Ftail) *(nz*Wo)^0.376*Q^0.122*S_vt^0.873*...
      ((100*t_c)/(cosd(Lambda_vt)))^(-0.49)*(AR_vt/(cosd(Lambda_vt)^2))^(0.357)* ...
      taper_vt^(0.039);


% Fuselage Weight
S_FUS = 257.6; % fuselage wetted area in ft^2 (SR-22)

W_FUS = 0.052 * S_FUS^(1.086) * (nz*Wo)^(0.177) * l_HT^(-0.051) * ...
   (l_FS/d_FS)^(-0.072) * Q^(0.241) + 11.9*(V_P*deltaP)^(0.271);


% Landing Gear Weight
W_LG = 0.025 * Wo;


% Avionics System Weight
W_AV = 50;


% Powerplant Weight
W_POWER = 1311.80;

% Fuel System Weight
W_FS = (M_h+M_o+M_ot+M_ht);


% Flight Control System Weight
W_CTRL = 0.053 * l_FS^(1.536) * b^(0.371) * (nz*Wo*10^-4)^(0.8);


% Electrical System Weight
W_EL = 12.57 * (W_FS+W_AV)^0.51;


% Air-Conditioning and Anti-icing
W_AC = 0.265 * Wo^0.52 * N_occ^0.68 * W_AV^(0.17) * M^(0.8);


% Furnishings
W_FURN = 0.0582 * Wo - 65;



% Total aircraft weight
W = Wing_WT + W_HT + W_VT + W_FUS + W_LG + W_CTRL +...
   W_EL + W_AC + W_FURN;


