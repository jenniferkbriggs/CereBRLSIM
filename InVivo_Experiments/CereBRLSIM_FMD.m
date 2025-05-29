
function [dy] = CereBRLSIM_FMD(t,y,paramvals,Pres,shear, etco2, eto2) %#codegen
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
%                                 satifies perfusion needs. If not present
%                                 the flow will set to 11.66
%   Pres        - Pressure [N/cm]
%   CBF         - From experiment

%   Learned parameters - 
%       Cmyo,   Cendo,   Cmeta  : Functionality of parameters
%       M                       : Metabolic state

%order: 
% 1:ropt, 2:cbf_init, 3:taumyo, 4:tauendo, 5:taumeta, 6:Cmyo, 7: Cendo, 8,
% Cmeta, 9: M, 10: baseline CO2, 11: eta, 12: alpha (controls gaussian), 13:
% beta_Length
%% ---------------------------------------------------------------------------  


%--- Pressure parameters ----- %
r0 = paramvals(1); %must be the same as y(2) initial
ropt = r0;
interpstep = 0.5; %interpolation step for interpolating Pressure driver

%--- Flow parameters --------- % 
Q0 = paramvals(2);


% --- K parameters ---- % (K defines the lenth-tension relationship of the vessel.)
K = @(r, ropt,alpha) 1/1.*exp(-(r-ropt).^2/alpha)+1e-8;
alphaK = r0./paramvals(12); %width of gaussian parameterized by the percent of the radius that can change


% --- Metabolic mechanism parameters --- %
%check inputs and make arrays:
if length(etco2)>0
    etco2 = Pres;
    etco2(:,2) = paramvals(10).*ones(length(Pres),1);
end

if length(eto2)>0
    po2 = eto2; %0 if we use simulated cbf, 1 if we use data.
    gamma = 70; %assume 10
    fo2 = 1;
    lambda = 1;%force starts at po2 \pm 5 mmhg above gamma
else %use metabolism;
    gamma = Q0/10;
    fo2 = 0;
    lambda = 10; %force starts at qn \pm 0.5 ml/s above gamma
end


% -- interpolate to current pressure and other drivers  ----- %
xq = [t-interpstep:1/10:t+interpstep]; %This is extrememly important and depends on the time step to be small
if t > Pres(end,1)+5 %for filtering
    xq = [Pres(end-1,1), Pres(end-1,1)+interpstep];
    disp('I do not know if this is right lets see if it ever gets called')
elseif t - paramvals(3) < Pres(1,1)
    xq = [min(t-interpstep,0):1/10:t+interpstep];
end

timeextra = 0.2;
Press_int = interp1(Pres(:,1), Pres(:,2),xq);
    while length(find(isnan(Press_int)))+3 > length(Press_int)
        xq = [t-paramvals(3)-timeextra:1/10:t+interpstep+timeextra];
        Press_int = interp1(Pres(:,1), Pres(:,2),(xq));
        timeextra = timeextra+0.2;
    end    
 Press_int = Press_int(~isnan(Press_int));
% %     
gradP = (diff(Press_int(end-2:end))/diff(xq(end-2:end))); %(phi, theta, t, tau, mi, ma);


% interpolate to shear data: 
if length(shear)>0
    xq = [t-interpstep:1/10:t+interpstep]; %This is extrememly important and depends on the time step to be small
    if t > shear(end,1)+5 %for filtering
        xq = [shear(end-1,1), shear(end-1,1)+interpstep];
        disp('I do not know if this is right lets see if it ever gets called')
    elseif t - paramvals(3) < shear(1,1)
        xq = [min(t-interpstep,0):1/10:t+interpstep];
    end
    
    timeextra = 0.2;
    CBF_init = interp1(shear(:,1), shear(:,2),xq);
        while length(find(isnan(CBF_init)))+3 > length(CBF_init)
            xq = [t-paramvals(3)-timeextra:1/10:t+interpstep+timeextra];
            CBF_init = interp1(shear(:,1), shear(:,2),(xq));
            timeextra = timeextra+0.2;
        end    
     CBF_init = CBF_init(~isnan(CBF_init));
    % %     
    Qdot = (diff(CBF_init(end-2:end))/diff(xq(end-2:end))); %(phi, theta, t, tau, mi, ma);
    noCBF = 0;
else
    noCBF = 1;
end

% ETCO2
if length(etco2>0)
e_init = interp1(etco2(:,1), etco2(:,2),xq);
    while length(find(isnan(e_init)))+3 > length(e_init)
        xq = [t-paramvals(3)-timeextra:1/10:t+interpstep+timeextra];
        e_init = interp1(etco2(:,1), etco2(:,2),(xq));
        timeextra = timeextra+0.2;
    end    
 e_init = e_init(~isnan(e_init));
% %     
e = e_init(1); %ml/s half of adequate CBF to initiate metabolism
else
    e = 40;
    etco2(1,2) = 40;
end

% PO2
if fo2 %if there are data for partial pressure of oxygen
    p_init = interp1(po2(:,1), po2(:,2),xq);
        while length(find(isnan(p_init)))+3 > length(p_init)
            xq = [t-paramvals(3)-timeextra:1/10:t+interpstep+timeextra];
            p_init = interp1(po2(:,1), po2(:,2),(xq));
            timeextra = timeextra+0.2;
        end    
     p_init = p_init(~isnan(p_init));
    % %     
    p = p_init(1); %ml/s half of adequate CBF to initiate metabolism
else
    p = y(1); %otherwise use flow
end


Kcurrent = K(y(2), ropt, alphaK);
% yyaxis left, plot(t,Kcurrent, 'o'), hold on, 
% yyaxis right, plot(t,y(2),'o')
% pause(0.2)
%---- temporal parameters ----- %
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

sendo = r0^3/Q0;

smyo = 1/r0.^4;
phi_myo = (smyo/taumyo)*Kcurrent*Cmyo*(P0*r0.^6-Press_int(1)*y(2).^6);

if exist('sh')
    h1_dot = 1/(2*tauendo) *((sh - shear(1,2))./shear(1,2));%y6
else
    h1_dot = 1/(2*tauendo) *((sendo*y(1)/y(2)^3 - sendo*Q0/r0^3));%y6
end
h2_dot = 1/(1*tauendo/2) *(y(6) - y(7)); 
h3_dot = 1/(2*tauendo/5) *(y(7) - y(8));
h4_dot = 1/(1*tauendo/10) *(y(8) - y(9));

fendo = (4*eta*r0/pi)*h4_dot;
phi_endo = Kcurrent*Cendo*fendo;


%metabolic
po2 = lambda*(p-gamma); 
M = paramvals(9);
xmeta = ((e-paramvals(10))/(paramvals(10)) + M) + 1/(1+exp(po2));
smeta = 1/100;
phi_meta = (smeta/taumeta)*(Kcurrent*Cmeta*xmeta);


%radius
r_dot = phi_myo+phi_endo+phi_meta;

if noCBF %if shear is not used as a driver in the model
    Qdot = alpha*gradP*(y(2))^4 + 4*alpha*(Press_int(1))*y(2)^3*r_dot;
end

% ---- Here's where the magic happens ------ %
dy = [Qdot;
     r_dot;
     phi_myo; %phi myo' : x(3) take m                                                                           ean Press_init to average over some hr.
     phi_endo; %phi endo' : x(4) 
     phi_meta; %phi meta' : x(5)
     %endo lag states
     h1_dot; %y(6)
     h2_dot; %y(7)
     h3_dot; %y(8)
     h4_dot]; %y(9) -- endo lag
end