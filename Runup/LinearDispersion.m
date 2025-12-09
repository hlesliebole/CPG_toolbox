function [L,C,Cg]=LinearDispersion(frequency,depth)

%% Linear Dispersion Equation Solver 

% Input: wave frequency (Hz), and water depth (m)

% Returns: wavelength L (m) , phase speed c (m/s), 
%          and group speed Cg (m/s) of a linear wave.

% Solves the linear dispersion relationship for surface gravity waves

%    (2pi/frequency)^2=g*(2pi/L)*tanh(depth*2pi/L)

%  for L, where g= gravitational acceleration (9.81m/s^2), 
%   using an algorithm derived by C.S. Wu.

% Wu, Chung-Shang & Thornton, Ed. (1986). Wave Numbers of Linear Progressive 
% Waves. J WATERW PORT COAST OC-ASCE. 112. 10.1061/(ASCE)0733-950X(1986)112:4(536). 

% This is one of many published approximate solutions. 

% For a detailed derivation and explanation of the linear dispersion
%  relationship see: http://falk.ucsd.edu/pdf/WavesLecture02_211A.pdf 

%% Apply the C.S. Wu algorithm to get radian wavenumber

a=4.0243*depth*frequency.^2;
yhat=a.*(1+1.26*exp((-1.84).*a));t=exp((-2).*yhat);
aka=a;aka(a >= 1)=a(a >= 1).*(1+2*t(a >= 1).*(1+t(a >= 1)));
aka(a < 1)=sqrt(a(a < 1)).*(1+(a(a < 1)./6).*(1+a(a < 1)./5));
k=abs(aka./depth); % radian wavenumber 2pi/L

% linear wave parameters

L=2*pi./k; % wavelength in meters
C=L.*frequency; % wave speed in m/s
% group velocity in m/s
sigma=2*depth.*k;
Cg=pi*(frequency./k).*(1+sigma./sinh(sigma));

end