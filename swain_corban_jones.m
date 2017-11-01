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

    function main
        cleanup;
        close all;
        figures = [fig1_protocol ...
            ];
        corban_figure_defaults;
        for fig = figures
            fh = makefigure(fig);
            savefig(fh);
        end
    end

%% Parameters
% Rate Constants and initial conditions - GET from Jones paper!
k1 = 2e-2;  % [1/nM/s] - Activation of V by Xa
k2 = 2e-2;  % [1/nM/s] - Activation of V by IIa
k3 = 1e-2;  % [1/nM/s] - Activation of VIII by Xa
k4 = 2e-3;  % [1/nM/s] - Activation of VIII by IIa (incorrect in paper)
k5 = 1e-2;  % [1/nM/s] - Conversion of mIIa to IIa by Va-Xa
k6 = 0.1;   % [1/nM/s] - on-rate for rapidly formed complexes
k7 = 1e-2;  % [1/nM/s] - on-rate for the VIIIa-IXa complex
k8 = 0.4;   % [1/nM/s] - on-rate for the Va-Xa complex
k9 = 5e-3;  % [1/s]   - off-rate for VIIIa-Ixa
k10 = 0.4;  % [1/s]   - off-rate for Va-Xa complex
k11 = 0.3;  % [1/s]   - Vmax for activation of IX by TF-VIIIa
k12 = 1.15; % [1/s]   - Vmax for activation of X by TF-VIIIa
k13 = 8.2;  % [1/s]   - Vmax for activation of X by VIIIa-IXa
k14 = 32;   % [1/s]   - Vmax for mIIa formation by Va-Xa
k15 = 1e-4; % [1/nM/s] - Activation of IX by Xa
k16 = 24;   % [1/s]   - off-rate for IX on TF-VIIa
k17 = 44;   % [1/s]   - off-rate for X on TF-VIIa complex
k18 = 1e-3; % [1/s]   - off-rate for X on VIIIa-IXa complex
k19 = 70;   % [1/s]   - off-rate for II on Va-Xa complex
k20 = 2e-2; % [1/s]   - constant for the slow degration of VIIIa-IXa
p_original = collect_params;

% Initial concentrations in [nM] (GET from Lawson et al. 1994)
TF_VIIa = 5e-3;     % species 1 
IX = 90;            % species 2
X = 170;            % species 3
V = 20;             % species 4
VIII = 0.7;         % species 5
II = 1.4e3;         % species 6 - prothrombin
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
y0_original = collect_initials;

    function p = collect_params
    % collect parameters for original model
        p = [k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,...
             k11,k12,k13,k14,k15,k16,k17,k18,k19,k20];
    end

    function y0 = collect_initials
    % Collect species initial concentrations for original model
        y0 = [TF_VIIa,IX,X,V,VIII,II,VIIIa_IXa,Va_Xa,IIa,Va_Xa_II,mIIa,...
            TF_VIIa_IX,TF_VIIa_X,VIIIa_IXa_X,IXa,Xa,Va,VIIIa];
    end

% ODE solver options
tol = 1e-9;
odeopts = odeset('RelTol',tol,'AbsTol',tol);
tspan = [0, 250];

%% Index Struct
% Struct for easier indexing of specific species.
ind.TF_VIIa = 1;
ind.IX =  2;
ind.X =  3;
ind.V =  4;
ind.VIII = 5;
ind.II = 6;
ind.VIIIa_IXa =  7;
ind.Va_Xa = 8;
ind.IIa = 9;
ind.Va_Xa_II = 10;
ind.mIIa = 11;
ind.TF_VIIa_IX = 12;
ind.TF_VIIa_X = 13;
ind.VIIIa_IXa_X = 14;
ind.IXa = 15;
ind.Xa = 16;
ind.Va = 17;
ind.VIIIa = 18;

%% Useful Calculations
thromb_form = @(y) y(:, ind.mIIa) * 1.2 + y(:, ind.IIa);

%% Figures
    function figspec = fig1_protocol
        p = p_original;
        y0 = y0_original;
        y0(ind.TF_VIIa) = 5e-3; % nM
        [t, y] = ode15s(@odefun, tspan, y0, odeopts, p);
        figspec.n = 1;
        figspec.title = '1A-Thrombin Timecourse, Model Validation';
        figspec.position = [3 640 606 324];
        figspec.x = t;
        figspec.y = thromb_form(y) / y0(ind.II) * 100;
        figspec.xlabel = 'Time (seconds)';
        figspec.ylabel = '% Thrombin Formation';
        figspec.xlim = [0 250];
        figspec.ylim = [0 120];
    end

    function figspec = fig2_protocol
        p = p_original;
        p([1, 2, 20]) = [100, 1e3, 1e5];
        y0 = y0_original;
        [t, y] = ode15s(@odefun,tspan,y0,odeopts,p);
    end

main;
end

function fighand = makefigure(fs)
fighand = setupfig(fs.n, fs.title, fs.position);
plot(fs.x, fs.y);
xlabel(fs.xlabel);
ylabel(fs.ylabel);
xlim(fs.xlim);
ylim(fs.ylim);
end

%% ODE Function
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

%% Corban Swain Utilities
function new_fig = setupfig(n,title,location)
% SETUPFIGURE Sets up a new figure.
%
% n: figure number
% title: figure title
% location: figure position, [left bottom width height]
%
new_fig = figure(n); clf; hold on;
box off;
new_fig.Name = title;
if nargin > 2
    new_fig.Position = location;
end
end

function savefig(fig,fig_name)
% SAVEFIGURE Saves the passed figure as a 300 dpi png.

if ~isdir([pwd filesep 'Figures'])
    mkdir 'Figures'
end
f = gobjects(1,1);
name = '';
switch nargin
    case 0
        f = gcf;
    case 1
        f = fig;
    case 2
        f = fig;
        name = fig_name;
end
if isempty(name)
    if isempty(f.Name)
        name = 'Untitled';
    else
        name = fig.Name;
    end
else
    if ~isempty(f.Name)
        name = [name, '-', f.Name];
    end
end
filename = ['Figures' filesep name];
print(f,filename,'-dpng','-r300');
end

function save_all_figures(trial_name)
% SAVEALLFIGURES saves all open figures as 300 dpi png files.

fprintf('%s - Saving all figures ...\n\n', datestr(now));
num = 1;
figs = findall(groot, 'Type', 'figure');
num_figs = length(figs);

for i = 1:num_figs
    f = figs(i);
    
    if isempty(f.Name)
        name = sprintf('Untitled%02d',num);
        num = num + 1;
    else
        name = f.Name;
    end
    
    if nargin == 1
        name = sprintf('%s - %s',trial_name,name);
    end
    
    fprintf('Saving Figure %d of %d, \"%s\" ... \n',...
        i, num_figs, name);
    savefig(f, name);
end
fprintf('Saving Done!\n\n')

end

function corban_figure_defaults
% CORBANFIGUREDEFAULTS Sets default values to make pretty figures.
fontSize = 15;
font = 'Helvetica';
set(groot, ...
    'defaultLineMarkerSize', 40,...
    'defaultLineLineWidth', 3, ...
    'defaultAxesFontSize', fontSize, ...
    'defaultAxesTitleFontWeight', 'normal', ...
    'defaultAxesFontName', font, ...
    'defaultAxesLabelFontSizeMultiplier', 1.1, ...
    'defaultAxesLineWidth', 2, ...
    'defaultFigureColor', [1 1 1], ...
    'defaultTextInterpreter', 'tex', ...
    'defaultTextFontSize',fontSize, ...
    'defaultTextFontName', font ...
    );
end

function cleanup
clc;
clear;
end