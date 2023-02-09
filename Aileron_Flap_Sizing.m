%% Aileron and Flap Sizing for GA Aircraft Wing

%% Aileron Sizing 

m_takeoff = 1134; % mass in kg
taper_ratio = 0.55; % Lambda
S_ht = 2.39; % Horizontal tail area in m^2
S_vt = 1.48; % Vertical tail area in m^2
V_stall = 28.29; % Stall speed in m/s
CLalpha = 5.4113; % lift-curve slope for NACA 644-421 (1/rad)
S = 11.6; % Wing area in m^2
AR = 10; % Aspect Ratio of wing
inner = 0.61; % inner location for the aileron
outer = 0.97; % outer location for the aileron
rho = 1.225; % sea level density in kg/m^3


% Step 1: Review rules/requirements for the aircraft you're designing

% Note that the maneuverability and roll control requirements for general 
% aviation aircraft can be found in FAR 23.147.


% Step 2: Choose roll control surface configuration


% Step 3: Refer to Figure 12.12 from the book to get the following:

tao = 0.41; % Control surface angle of attack effectiveness parameter
%C_A = 0.2*C; % Chord length of ailerons relative to the wing chord

% Step 4: Determine the Class, Flight Phase, and Level of Acceptance required

%Class = I;
%Flight_Phase = B;
%Lvl_of_Accpt = 1;


% Step 5: Roll control handling qualities design requirements:
Phi_desired = 45 * (pi/180); % required bank angle in radians
t_bank45_req = 1.7; % Time to achieve bank angle of 45 deg in sec


% Step 6: Assuming flap location extends to 60% of wingspan, choose span 
% of aileron.

%bai = 0.7*b; % interior edge location along wing of the aileron
%bao = 0.95*b; % exterior edge location along wing of the aileron


% Step 9: Find aileron rolling moment coefficient derivative

b = (AR*S)^(1/2); % span of wing
Cbar = b/AR; % Average chord length along wing
Cr = Cbar*(3/2)*((1+taper_ratio) / (1 + taper_ratio + taper_ratio^2)); % wing root chord 
I_xx = (1/12)*m_takeoff*(b/2)^2; % Rolling moment of inertia in kg*m^2


y_i = inner*(b/2); % interior edge location along wing of the aileron
y_o = outer*(b/2); % exterior edge location along wing of the aileron

y_D = 0.4*(b/2); % Drag moment arm (estimated as 40% of wing span)

CldeltaA = ((2*CLalpha*tao*Cr)/(S*b)) * ( (y_o^2 / 2 + (2/3)*((taper_ratio-1) / b)*y_o^3) - (y_i^2 / 2 + (2/3)*((taper_ratio-1) / b)*y_i^3) );

% Step 10: Assume a maximum deflection angle
Delta_A = 15; % Max degrees +/-

% Step 11: Calculate aircraft rolling moment coefficient when ailerons are at max
% deflection 
Cl = CldeltaA*Delta_A * (pi/180);

V_approach = 1.3*V_stall; % Approach speed in m/s
L_A = (1/2)*rho*V_approach^2*S*Cl*b; % rolling moment coefficient at max deflection angle

% Step 13: Calculate the steady state roll rate
CDR = 0.9; % wing-horizontal tail-vertical tail rolling drag coefficient 
P_ss = sqrt((2*L_A) / (rho * (S+S_vt+S_ht) * CDR * y_D^3)); % rad/sec

% Step 14: Calculate the bank angle where the aircraft achieves the steady
% state roll rate:

Phi_1 = I_xx / (rho*y_D^3 * (S+S_vt+S_ht) * CDR) * log(P_ss^2); % rad

% Step 15: Calculate the aircraft rate of roll produced by the aileron
% rolling moment until aircraft reaches the steady state roll rate

P_dot = (P_ss)^2/(2*Phi_1);

% Steps 16 and 17: Compare bank angle from step 5 to bank angle in step 14
% and calculate time required to bank 45 degrees. Compare this to the
% required time for GA aircraft outlined in step 5.

time = sqrt((2*Phi_desired) / (P_dot)); % required time in sec

% Interpolate to find the chord length at each end of the aileron
C_tip = (taper_ratio)*(Cr);
x1 = 0;
y1 = Cr;
x2 = b/2;
y2 = C_tip;

C_inner = (y2-y1)/(x2-x1)*(y_i-x1)+y1; 
C_outer = (y2-y1)/(x2-x1)*(y_o-x1)+y1; 

% Dimensions of the ailerons
Aileron_inner_C = 0.2*C_inner;
Aileron_outer_C = 0.2*C_outer;
b_Aileron = (outer*(b/2) - inner*(b/2));

fprintf('The inner aileron chord-length is %1.3f meters.\n\n', Aileron_inner_C);
fprintf('The outer aileron chord-length is %1.3f meters.\n\n', Aileron_outer_C);
fprintf('The span of each aileron is %1.3f meters.\n\n', b_Aileron);

%% Flap Sizing

m_takeoff = 1134; % mass in kg
taper_ratio = 0.55; % Lambda
S_ht = 2.39; % Horizontal tail area in m^2
S_vt = 1.48; % Vertical tail area in m^2
V_stall = 28.29; % Stall speed in m/s
Cr = 1.3517; % wing root chord
CLalpha = 5.4113; % lift-curve slope for NACA 644-421 (1/rad)
S = 11.6; % Wing area in m^2
AR = 10; % Aspect Ratio of wing
b = (AR*S)^(1/2); % span of wing
inner = 0; % inner location for the flap
outer = 0.6; % outer location for the flap
rho = 1.225; % sea level density in kg/m^3
g = 9.81; % gravitational constant in m/s^2


V_TO = 1.2*V_stall; % average takeoff speed

CL_TO = 0.85 * (2*m_takeoff*g) / (rho*(V_TO)^2*S); % wing lift coeffient required at takeoff
% NOTE: the 0.85 factor is because during take-off, the aircraft has a
% take-off angle (say about 10 deg). Thus about 15% of the lift is maintained by the
% vertical component (sin (10)) of the engine thrust.

% CLwing = ; % Total lift coefficient generated by wing at takeoff (Lifting Line Theory w/o flap)
% CL_TO_actual = ; % Total lift coefficient generated (w/ flaps) for takeoff


% Interpolate to find the chord length at each end of the aileron
C_tip = (taper_ratio)*(Cr);
x1 = 0;
y1 = Cr;
x2 = b/2;
y2 = C_tip;

y_i = inner*(b/2); % interior edge location along wing of flap
y_o = outer*(b/2); % exterior edge location along wing of flap

C_inner = (y2-y1)/(x2-x1)*(y_i-x1)+y1; 
C_outer = (y2-y1)/(x2-x1)*(y_o-x1)+y1; 

% Dimensions of the ailerons
Flap_inner_C = 0.2*C_inner;
Flap_outer_C = 0.2*C_outer;
b_Flap = (outer*(b/2) - inner*(b/2));

fprintf('\n\nThe inner flap chord-length is %1.3f meters.\n\n', Flap_inner_C);
fprintf('The outer flap chord-length is %1.3f meters.\n\n', Flap_outer_C);
fprintf('The span of each flap is %1.3f meters.\n\n', b_Flap);




