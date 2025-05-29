
function [dy] = CereBRLSIM_dCA(t,y,paramvals,Pres, etco2, pao2, pbto2) %#codegen
%% -----------------------------------------------------------------------------
% Cerebral Blood Flow with temporaly informed cerebral vascular regulation paramters 
% Model is original but was motivated by Ursino and Lodi 1998
% Jennifer Briggs 2022
%% -----------------------------------------------------------------------------


%% ---------------------------------------------------------------------------  
% INPUTS:
%   t           - time comes from data
%   y           - system of 6 ODEs
%   paramvals   - 
%       taumyo, tauendo, taumeta: time lag (s) for three mechanisms
%       ropt                    : radius for which the force is maximal (affects K)
%       cbf_init                : Value of inital blood flow. This is
%                                 assumed to be `normal flow' which
%                                 satifies perfusion needs. 
%   Pres        - Pressure [N/cm]

%   OPTIONAL INPUTS: 2xT, where first column is time and second column is
%   data
%   etco2       - End-tital CO2 From ventillator [mmHg]
%   pao2        - Partial pressure of arterial oxygen - from data [mmHg]
%   pbto2       - Partial pressure of brain tissue oxygen - from data
%   [mmHg]


%   Learned parameters - 
%       Cmyo,   Cendo,   Cmeta  : Functionality of parameters
%       qn                      : Perfusion demand of brain/onset of
%                               ischemia
%       bco2                    : Sensativity of the metabolic mechanism.

%order: 
% 1:r0, 2:cbf_init, 3:taumyo, 4:tauendo, 5:taumeta, 6:Cmyo, 7: Cendo, 8,
% Cmeta, 9: M, 10: etco2_baseline, 11: eta, 12: alpha (controls gaussian), 13:
% beta_Length, 14: ropt, 15: Qopt, 16: Popt
%% ---------------------------------------------------------------------------  


%--- Pressure parameters ----- %
%if ropt is manually set to be different than r0
if length(paramvals) > 13 & paramvals(14) ~= 0
    ropt = paramvals(14);
else
    ropt = paramvals(1);
end

%per = 1./mean(diff(Pres(1,:)));
interpstep = 0.5; %interpolation step for interpolating Pressure driver

if length(paramvals) > 14 & paramvals(15) ~= 0
    Q0 = paramvals(15);
else
    Q0 = paramvals(2);
end

P0 = mean(Pres(1:5,2));


per = 1./mean(diff(Pres(1,:)));
interpstep = 0.5; %interpolation step for interpolating Pressure driver
r0 = paramvals(1);



% --- K parameters ---- % (K defines the lenth-tension relationship of the
% vessel. i.e. how much force the vessel)
K = @(r, ropt,alpha) 1/1.*exp(-(r-ropt).^2/alpha)+1e-8;
alphaK = r0./paramvals(12); %width of gaussian parameterized by the percent of the radius that can change


% --- Metabolic mechanism parameters --- %
%check inputs and make arrays:
if ~exist('etco2')
    etco2 = Pres;
    etco2(:,2) = 40.*ones(length(Pres),1);
end

if exist('pao2')
    po2 = pao2; %0 if we use simulated cbf, 1 if we use data.
    gamma = 50;
    fo2 = 1;
    lambda = 1;%force starts at po2 \pm 5 mmhg above gamma
elseif exist('pbto2')
    po2 = pbto2;
    gamma = 10;
    fo2 = 1;
    lambda = 5; %force starts at po2 \pm 1 mmhg above gamma
else %use metabolism; 
    gamma = Q0/10;
    fo2 = 0;
    lambda = 10; %force starts at qn \pm 0.5 ml/s above gamma
end




% try
% -- solve for P and K  ----- %
xq = [t-interpstep:1/10:t+interpstep]; %This is extrememly important and depends on the time step to be small
if t > Pres(end,1)+5 %for filtering
    xq = [Pres(end-1,1), Pres(end-1,1)+interpstep];
    disp('I do not know if this is right lets see if it ever gets called')
elseif t - paramvals(3) < Pres(1,1)
    xq = [t-interpstep:1/10:t+interpstep];
end

    Press_int = interp1(Pres(:,1), Pres(:,2),xq);
  
    Press_int = Press_int(~isnan(Press_int));
    gradP = (diff(Press_int(1:3))/diff(xq(1:3))); %(phi, theta, t, tau, mi, ma);


% ETCO2
e_init = interp1(etco2(:,1), etco2(:,2),xq); 
e_init = e_init(~isnan(e_init));
% %     
e = mean(e_init(1:3)); %ml/s half of adequate CBF to initiate metabolism


% PO2
if fo2 %if there are data for partial pressure of oxygen
    p_init = interp1(po2(:,1), po2(:,2),xq);
     p_init = p_init(~isnan(p_init));
    % %     
    p = mean(p_init(1:3)); %ml/s half of adequate CBF to initiate metabolism
else
    p = y(1); %otherwise use flow
end


Kcurrent = K(y(2), ropt, alphaK);


%K may need to be positive - check signs

% ---- temporal parameters ----- %
taumyo = paramvals(3);
tauendo = paramvals(4);
taumeta = paramvals(5);


% ---- Inferred parameters ----- %
Cmyo = paramvals(6);
Cendo = paramvals(7);
Cmeta = paramvals(8);
M = paramvals(9);

% --------- Other Parameters ---------- %
eta = paramvals(11);

alpha = pi/(8*paramvals(11)*paramvals(1).^4./paramvals(13)); %set alpha
P0 = Pres(1,2);%original pressure 

%maybe(?) find a scaling value for xQ/r0^3 = 1
x = r0^3/Q0;

%define states:
%myogenic
smyo = 1/r0.^4;
phi_myo = (smyo/taumyo)*Kcurrent*Cmyo*(P0*r0.^6-Press_int(1)*y(2).^6);

%endothelial
h1_dot = 1/(1*tauendo/2) *((x*y(1)/y(2)^3 - x*Q0/r0^3));%y6
h2_dot = 1/(1*tauendo/5) *(y(6) - y(7)); 
h3_dot = 1/(2*tauendo/10) *(y(7) - y(8));
h4_dot = 1/(1*tauendo/10) *(y(8) - y(9));
fendo = (4*eta*r0/pi)*h4_dot;
phi_endo = Kcurrent*Cendo*fendo;
 

%metabolic
po2 = lambda*(p-gamma);
xmeta = ((e-paramvals(10))/paramvals(10) + M+ 1/(1+exp(po2)));
smeta = 1/100;
phi_meta = (smeta/(taumeta))*(Kcurrent*Cmeta*xmeta);

% %radius
r_dot = phi_myo+phi_endo+phi_meta;


%plot(Press_int(1)*75, Kcurrent, 'ko'), hold on, drawnow
% ---- Here's where the magic happens ------ %
dy = [alpha*gradP*(y(2))^4 + 4*alpha*(Press_int(1))*y(2)^3*r_dot;
     r_dot;
     phi_myo; %phi myo' : x(3) take m                                                                           
     phi_endo; %phi endo' : x(4) 
     phi_meta; %phi meta' : x(5)
     h1_dot; %y(6)
     h2_dot; %y(7)
     h3_dot;
     h4_dot]; %y(8) -- endo lag
    %The mutiply by round(x*1000)/1000 allows us to convert to Mex it is
    %equivalent to rounding to the 4th digit

end