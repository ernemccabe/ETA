% Calculates Elliptical Lift Distribution For Our Aircraft

% Define Initial Conditions
AR=10; % aspect ratio
b=11.6586; % span [m] span is 38 ft 4 inches for SR22
WL=1197.01; % wing loading [N/m^2]
V=77.16667; % flying speed  at sea level [m/s]
rho=1.226; 
q=.5*rho*V^2; % [N/m^2]
m_0=2*pi;`

% projected area and total lift of aircraft are:
S=b^2/AR; % [m^2]
L=WL*S; % [N]

% using Eqns 6.21 and 6.23 to find wing lift and drag coefficients
Cl=L/(q*S);
Cd=Cl^2/(pi*AR);

% induced drag
D=L*Cd/Cl;

% Percent of drag compared to total lift
Percent=D/L*100;

% Additional power required to compensate for induced drag 
P=D*V;

% Induced angle of attack (Eqn 6.23) and the constant down wash (Eqn 6.4)
alpha_i=-Cd/Cl;
w=alpha_i*V;

% effective angle of attach and absolute angle of attack are constant along
% span and computed iwth Eqn 6.10 and 6.9
alpha_0=Cl/m_0;
alpha_a=alpha_0-alpha_i;

% lift curve slope defined in Eqn 6.11 and computed with Eqn 6.12
m=m_0/(1-alpha_i/alpha_0);