function out = closure_temp(d,cr)
%%Input definitions
% d = the characteristic diffusion distance (radius of sphere, half-width
% of cylinder, etc.)(um)
% cr = cooling rate (oC/Myr)

%%Constants
A = 27; %Dodson (1973) geometric factor; use 55 for sphere, 27 for cylinder, or 9 for sheet
R = .001986; %(kcal/(K*mol)); Gas constant

% Diffusivity values for Rutile (Cherniak, 1993)
% Ea_rtl = 89.2; %(kcal/mol); activation E for Pb in rt
% Do_rtl = 1.7e9; %((m^2)/Ma)diffusivity constant
% Diffusivity values for Apatite (Cherniak et al., 1991)
Ea = 54.6; %(kcal/mol); activation E for Pb in ttn
Do = 4.01e5; %((m^2)/Ma)diffusivity constant

%convert um to m
d = d*0.000001;

syms Tc
PbTc = solve(Tc==Ea/(R*log((A*(((Tc^2)*R)/(Ea*cr))*Do)/(d^2))),Tc);

%Output is degrees Celcius
out = PbTc-273;
end