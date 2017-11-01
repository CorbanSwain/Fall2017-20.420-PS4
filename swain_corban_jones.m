function swain_corban_jones
% Implementation of Tissue Factor Pathway to Thrombin model for course
% 20.420 - October, 2015
%
% Reference: Jones, KC and Mann, KG. A model for the tissue factor pathway
% to thrombin. J Biol Chem 269:37, 1994.
%
% Initial conditions reference: Lawson, JH, et al. A model for the tissue
% factor pathway to thrombin. J Biol Chem 269:37, 1994.
%
% HW Completed by Corban Swain
% November 2017

% FIXME - Remove Close All before submission
clc; clear; % close all;
% Rate Constants and initial conditions - GET from Jones paper!
k1 = 2e7;   % [1/M/s] - Activation of V by Xa
k2 = 2e7;   % [1/M/s] - Activation of V by IIa
k3 = 1e7;   % [1/M/s] - Activation of VIII by Xa
k4 = 2e6;   % [1/M/s] - Activation of VIII by IIa (incorrect in paper)
k5 = 1e7;   % [1/M/s] - Conversion of mIIa to IIa by Va-Xa
k6 = 1e8;   % [1/M/s] - on-rate for rapidly formed complexes
k7 = 1e7;   % [1/M/s] - on-rate for the VIIIa-IXa complex
k8 = 4e8;   % [1/M/s] - on-rate for the Va-Xa complex
k9 = 5e-3;  % [1/s]   - off-rate for VIIIa-Ixa
k10 = 0.4;  % [1/s]   - off-rate for Va-Xa complex
k11 = 0.3;  % [1/s]   - Vmax for activation of IX by TF-VIIIa
k12 = 1.15; % [1/s]   - Vmax for activation of X by TF-VIIIa
k13 = 8.2;  % [1/s]   - Vmax for activation of X by VIIIa-IXa
k14 = 32;   % [1/s]   - Vmax for mIIa formation by Va-Xa
k15 = 1e5;  % [1/M/s] - Activation of IX by Xa
k16 = 24;   % [1/s]   - off-rate for IX on TF-VIIa
k17 = 44;   % [1/s]   - off-rate for X on TF-VIIa complex
k18 = 1e-3; % [1/s]   - off-rate for X on VIIIa-IXa complex
k19 = 70;   % [1/s]   - off-rate for II on Va-Xa complex
k20 = 2e-2; % [1/s]   - constant for the slow degration of VIIIa-IXa
    
% Initial concentrations in [M] (GET from Lawson et al. 1994)
TF_VIIa = 1;        % species 1 
IX = 90e-9;         % species 2
X = 170e-9;         % species 3
V = 20e-9;          % species 4
VIII = 0.7e-9;      % species 5
II = 1.4e-6;        % species 6 - prothrombin
VIIIa_IXa = 0;      % species 7
Va_Xa = 0;          % species 8
IIa = 0;            % species 9 - alpha-thrombin
Va_Xa_II = 0;       % species 10
mIIa = 0;           % species 11 - meizothrombin 
TF_VIIa_IX = 0;     % species 12
TF_VIIa_X = 0;      % species 13
VIIIa_IXa_X = 0;    % species 14
IXa = 0;            % species 15
Xa = 0;             % species 16
Va = 0;             % species 17
VIIIa = 0;          % species 18

% ODE solver options
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
tspan = [0 20];

%% Fig 1

% collect parameters for original model
p = [k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,...
     k11,k12,k13,k14,k15,k16,k17,k18,k19,k20];

% Collect species initial concentrations for original model
y0 = [TF_VIIa,IX,X,V,VIII,II,VIIIa_IXa,Va_Xa,IIa,Va_Xa_II,mIIa,...
      TF_VIIa_IX,TF_VIIa_X,VIIIa_IXa_X,IXa,Xa,Va,VIIIa];
  
[t1,y1] = ode15s(@odefun,tspan,y0,options,p);

%% Fig 2

k1 = 100;
k2 = 1000;
k20 = 10000;

% collect parameters for original model
p = [k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,...
     k11,k12,k13,k14,k15,k16,k17,k18,k19,k20];
  
[t2,y2] = ode15s(@odefun,tspan,y0,options,p);

% Call plotting of figure 1
FIG1(t1,y1);
FIG2(t2,y2);

end

function FIG1(t,y)
size(y)

figure(1); clf;
plot(t,y(:,1),t,y(:,2),t,y(:,6));
xlabel('Time (min)')
ylabel('Species (M)')

end

function FIG2(t,y)

figure(2); clf;
plot(t,y); 
xlabel('Time (min)')
ylabel('Species (M)')

end

% ODE FUNCTION
function ydot = odefun(~,y,p)

% Collect param values in cell array and redefine params with names
paramsCell = num2cell(p);
[k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,...
 k11,k12,k13,k14,k15,k16,k17,k18,k19,k20]=paramsCell{:};

% Collect y-vals in cell array and redefine y-vals with names
yCell = num2cell(y);
[TF_VIIa,IX,X,V,VIII,II,VIIIa_IXa,Va_Xa,IIa,Va_Xa_II,mIIa,...
 TF_VIIa_IX,TF_VIIa_X,VIIIa_IXa_X,IXa,Xa,Va,VIIIa] = yCell{:};

%%% ODEs %%%

% TF-VIIa - species 1
dTF_VIIa = k11*TF_VIIa_IX - k6*TF_VIIa*IX + ...
           k16*TF_VIIa_IX + k12*TF_VIIa_X -k6*TF_VIIa*X + k17*TF_VIIa_X;

% IX - species 2
dIX = k16*TF_VIIa_IX - k6*TF_VIIa*IX - k15*IX*Xa - k15*IX*Va_Xa;

% X - species 3
dX = k17*TF_VIIa_X - k6*TF_VIIa*X - k6*VIIIa_IXa*X + k18*VIIIa_IXa_X;    

% V - species 4
dV = -k1*V*Xa - k2*V*IIa - k2*V*mIIa;

% VIII - species 5
dVIII = -k3*VIII*Xa - k4*VIII*IIa - k4*VIII*mIIa;

% II - species 6
dII = k19*Va_Xa_II - k6*Va_Xa*II;

% VIIIa-IXa (first-order deg) - species 7 
dVIIIa_IXa = k7*VIIIa*IXa - k9*VIIIa_IXa - k6*VIIIa_IXa*X + ...
             k18*VIIIa_IXa_X + k13*VIIIa_IXa_X - k20*VIIIa_IXa;

% Va-Xa (incorrect in original paper) - species 8
dVa_Xa = k8*Xa*Va - k10*Va_Xa + k19*Va_Xa_II - k6*Va_Xa*II + k14*Va_Xa_II;

% IIa - species 9
dIIa = k5*Va_Xa*mIIa;

% Va-Xa-II - species 10
dVa_Xa_II = k6*Va_Xa*II - k19*Va_Xa_II - k14*Va_Xa_II;

% mIIa - species 11
dmIIa = k14*Va_Xa_II - k5*Va_Xa*mIIa;

% TF-VIIa-IX - species 12
dTF_VIIa_IX = k6*TF_VIIa*IX - k16*TF_VIIa_IX - k11*TF_VIIa_IX;

% TF-VIIa-X - species 13
dTF_VIIa_X = k6*TF_VIIa*X - k17*TF_VIIa_X - k12*TF_VIIa_X;

% VIIIa_IXa_X - species 14
dVIIIa_IXa_X = k6*VIIIa_IXa*X - k18*VIIIa_IXa_X - k13*VIIIa_IXa_X;

% IXa - species 15
dIXa = k9*VIIIa_IXa - k7*VIIIa*IXa + k11*TF_VIIa_IX + k15*(IX*Xa + IX*Va_Xa);

% Xa (incorrect in original paper: k8 should be used) - species 16
dXa = k10*Va_Xa - k8*Xa*Va + k12*TF_VIIa_X + k13*VIIIa_IXa_X;

% Va (incorrect in original paper: k8 should be used) - species 17
dVa = k10*Va_Xa - k8*Xa*Va + k1*V*Xa + k2*V*IIa + k2*V*mIIa;

% VIIIa - species 18
dVIIIa = k9*VIIIa_IXa - k7*VIIIa*IXa + k3*VIII*Xa + k4*(VIII*IIa + VIII*mIIa);

% Collect all ODEs for output
ydot = [dTF_VIIa;dIX;dX;dV;dVIII;dII;dVIIIa_IXa;dVa_Xa;dIIa;dVa_Xa_II;
        dmIIa;dTF_VIIa_IX;dTF_VIIa_X;dVIIIa_IXa_X;dIXa;dXa;dVa;dVIIIa];
end