% Calculates the center of mass for component placement in the aircraft

% Mass of hygrogen tanks
M_b=52.895; % [lb] Mass of large hydrogen tank
V_b=9.239; % [ft^3]Volume of large hydrogen tank (NOTE: taken from sheets, not same as one given in discord)
rho=M_b/V_b; % [lb/ft^3] density of material assuming uniform denisty
V_s=pi*(9.84/12)^2*(26.26/12); % [ft^3] v=pi*r^2*h:volume of small tank
M_s=rho*V_s; % [lb] mass of small tank assuming uniform density

% Masses from Powerplant Spreadsheet
M1=800.052; % [lb] Fuel Cell Stack
M2=22.04; % [lb] Weight of O2 tank 
M2_2=22.04+11.02; %[lb] Weight of O2 tank and O2
M3=116.812; % [lb] Weight of Battery
M4=2*M_s; % [lb] Weight of both H2 tanks beacuse they are in the same position
M4_2=2*(M_s+32.619); % [lb] Weight of H2 tanks and H2
M5=198.36; % [lb] Wieght of motor assembly
M6=78; % [lb] Weight of propeller
M7=51.3861; % [lb] Weight of nose section
M8=97.3858; % [lb] Weight of cabin section
M9=133.4139; % [lb] Weight of tail section
M10=623.71; % [lb] Weight of wings
M11=38.5+61; % [lb] Weight of Tail

% Distances from quarter chord with quarter chord being x=0
x1=8.43/12; % [ft] Distance from fuel cell stack to 1/4 chord
x2=-75.26/12; % [ft] Distance from O2 tank to 1/4 chord
x3=-73.91/12; % [ft] Distance from battery to 1/4 chord
x4=-54.22/12; % [ft] Distance from H2 tanks to 1/4 chord
x5=-71.54/12; % [ft] Distance from motor to 1/4 chord
x6=-110.28/12; % [ft] Distance from nose to 1/4 chord
x7=(-22.17-80.37)/12; %[ft] Distance from nose center of mass to 1/4 chord
x8=(42.58-80.37)/12; % [ft] Distance from cabin center of mass to 1/4 chord
x9=(142.72-80.37)/12; % [ft] Distance from tail center of mass to 1/4 chord
x10=0; % [ft] Distance from wing to 1/4 chord
x11=214.03/12; % [ft] Distance from Tail to 1/4 chord

% Calculate the x location of the center of maxx
x_CM_empty=(M1*x1+M2*x2+M3*x3+M4*x4+M5*x5+M6*x6+M7*x7+M8*x8+M9*x9+M10*x10+M11*x11)/(M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+M11); % [ft]
x_CM_full=(M1*x1+M2_2*x2+M3*x3+M4_2*x4+M5*x5+M6*x6+M7*x7+M8*x8+M9*x9+M10*x10+M11*x11)/(M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+M11); % [ft]

