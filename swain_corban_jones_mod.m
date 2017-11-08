function swain_corban_jones_mod
% Modification of implementation of Tissue Factor Pathway to Thrombin
% model for course 20.420 - November, 2017
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
        cleanup; % close all;
        fprintf('Beginning Script ...\n');
        figs_to_plot = 1:3;
        corban_figure_defaults;
        for i = figs_to_plot
            fs = figures.(sprintf('f%d', i))();
            fprintf('Drawing Figure %2d\n', i);
            fh = makefigure(fs); figure(fh);
%             savefig(fh);
        end
        fprintf('Done!\n');
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
k9 = 5e-3;  % [1/s]   - off-rate for VIIIa-IXa
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
km1 = 8.1e-3;
km1_2 = km1*0.1;
km2 = 1e-2;
km3 = 0.5e-3;
km4 = km2;
km5 = km3;
km6 = 1e-2;
km7 = 1e-2;
km8 = 8e-4;
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
I = 5e-3;           % upper limit on factor VIIIa_IXa
S = 332;
PC = 65;
APC = 0;
APC_S = 0;
APC_S_V = 0;
y0_original = collect_initials;

    function p = collect_params
    % collect parameters for original model
        p = [k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,...
             k11,k12,k13,k14,k15,k16,k17,k18,k19,k20, ...
             km1,km1_2,km2,km3,km4,km5,km6,km7,km8];
    end

    function y0 = collect_initials
    % Collect species initial concentrations for original model
        y0 = [TF_VIIa,IX,X,V,VIII,II,VIIIa_IXa,Va_Xa,IIa,Va_Xa_II,mIIa,...
            TF_VIIa_IX,TF_VIIa_X,VIIIa_IXa_X,IXa,Xa,Va,VIIIa,I,S,PC,...
            APC,APC_S,APC_S_V];
    end

% ODE solver options
tol = 1e-9;
odeopts = odeset('RelTol',tol,'AbsTol',tol,...
    'NonNegative',1:length(y0_original));
tspan = [0, 1e4];

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
ind.I = 19;
ind.S = 20;
ind.PC = 21;
ind.APC = 21;
ind.APC_S = 22;
ind.APC_S_V = 23;

%% Useful Calculations
thromb_activity = @(y) sum(y(:,[ind.mIIa, ind.IIa]) .* [1.2, 1], 2);
thromb_percent = @(y, II_0) thromb_activity(y) / II_0 * 100;
odesim = @(y0, p) ode15s(@odefun, tspan, y0, odeopts, p);

    function ps_out = plotdefaults(ps_in)
        ps_in.xlim = tspan / 60;
        ps_in.xlabel = 'Time (min)';
        ps_out = ps_in;
    end

    function [t, y, y0, p] = sim_from_maps(initial_maps, param_maps)
        % Runs a series of simulations using various changes as specified
        % by matrices having shape (n, 2). The first column is the index
        % of the species/param to be changed and the second column 
        % indicates the new values.
        n = length(initial_maps);
        if n ~= length(param_maps)
            error('Initial and param maps must be of the same size');
        end
        t = cell(n, 1); y = t; y0 = t; p = t;
        for i = 1:n
            y0{i} = y0_original;
            p{i} = p_original;
            im = initial_maps{i};
            if size(im,1) > 0
                y0{i}(im(:,1)) = im(:, 2);
            end
            pm = param_maps{i};
            if size(pm,1) > 0 
                p{i}(pm(:,1)) = pm(:,2);
            end
            [t{i}, y{i}] = odesim(y0{i}, p{i});
        end
    end

%% Figures
figures = struct;

figures.f1 = @fig1;
    function fs = fig1
        % Comparing effects of adding anticoagulatory factors
        fignum = 1;
        fprintf('Running Figure %2d\n', fignum);
        ntrials = 2;
        initial_map = cell(1, ntrials);
        param_map = initial_map;
        param_map{1} = [(21:29)', zeros(9, 1)];
        
        [t, y, y0] = sim_from_maps(initial_map, param_map);
        for i = 1:ntrials
            ps.x{i} = t{i};
            ps.y{i} = thromb_percent(y{i}, y0{i}(ind.II));
        end
        ps = plotdefaults(ps);
        ps.ylabel = '% Thrombin Formation';
        ps.ylim = [0 120];
        ps.xlim = [0, 250];
        ps.legend = {'Initial Model'; 'Anticoagulation Model'};
        ps.legend_loc = 'east';

        fs.n = fignum;
        fs.title = sprintf('MOD%d - Anti Coagulation Model', ...
            fignum);
        fs.position = [1000 548 929 407];
        fs.plots = {ps};
        fs.sub = [1, 1];
    end

figures.f2 = @fig2;
    function fs = fig2
        % Comparing Different contraceptive effects
        fignum = 2;
        fprintf('Running Figure %2d\n', fignum);
        ntrials = 4;
        initial_map = cell(1, ntrials);
        param_map = initial_map;
        species_incr = ["V";"VIII";"IX";"X";"PC"];
        n1 = length(species_incr);
        initial_map{2} = zeros(n1, 2);
        for i = 1:n1
            index = ind.(species_incr(i));
            initial_map{2}(i, :) = [index, y0_original(index) * 1.3];
        end
        initial_map{3} = [initial_map{2}; ...
            ind.S, y0_original(ind.S) * 0.8];
        initial_map{4} = [initial_map{3}; ...
            ind.TF_VIIa, y0_original(ind.TF_VIIa) * 1.30];
        
        [t, y, y0] = sim_from_maps(initial_map, param_map);
        for i = 1:ntrials
            ps.x{i} = t{i} / 60;
            ps.y{i} = thromb_percent(y{i}, y0{i}(ind.II));
        end
        ps = plotdefaults(ps);
        ps.ylabel = '% Thrombin Formation';
        ps.ylim = 'auto';
        ps.legend = {'None'; 'Progestrin'; 'Levono./Noreth.'; ...
            'Deso./Gesto.'};
        ps.legend_title = 'Contrtaceptive';
        ps.legend_loc = 'southeast';

        fs.n = fignum;
        fs.title = sprintf('MOD%d - Effects of Contraceptives', ...
            fignum);
        fs.position = [3 548 929 407];
        fs.plots = {ps};
        fs.sub = [1, 1];
    end

figures.f3 = @fig3;
    function fs = fig3
        fignum = 3;
        fprintf('Running Figure %2d\n', fignum);
        % Varying Extrinsic Pathway Activation
        TF_VIIa_initials = [0.1 5]; % nM
        ninits = length(TF_VIIa_initials) + 1;
        ntrials = ninits * 4;
        initial_map = cell(1, ntrials);
        param_map = initial_map;
        for i = 1:(ninits - 1)
            initial_map{i + 1}  = [ind.TF_VIIa, TF_VIIa_initials(i)];
        end
        species_incr = ["V";"VIII";"IX";"X";"PC"];
        n1 = length(species_incr);
        for i = (1:ninits) + ninits
            initial_map{i} = zeros(n1, 2);
            for j = 1:n1
                index = ind.(species_incr(j));
                initial_map{i}(j, :) = [index, y0_original(index) * 1.3];
            end
            initial_map{i} = [initial_map{i - ninits}; initial_map{i}];
        end
        for i = (1:ninits) + (ninits * 2)
            initial_map{i} = [initial_map{i - ninits}; ...
                ind.S, y0_original(ind.S) * 0.8];
        end
        for i = (1:ninits) + (ninits * 3)
            initial_map{i} = [initial_map{i - ninits}; ...
                ind.TF_VIIa, y0_original(ind.TF_VIIa) * 1.30];
        end
        
        
        [t, y] = sim_from_maps(initial_map, param_map);
        
        ps = struct;
        ps = plotdefaults(ps);
        ps.color_cycle = ninits;
        ps.spec_cycle = {'-';'--';':';'-.'}; 
        ps.x = t;
        for i = 1:ntrials
            ps.y{i} = thromb_activity(y{i}) ./ 1e3; 
        end
        ps.legend = {'Initial Model (5 pM)', ...
            '100 pM', '1 nM'};
        ps.legend_loc = 'southeast';
        ps.legend_title = '[TF-VIIa]_{0}';
        ps.ylabel = 'Thrombin Formation (\muM)';
        ps.ylim = [0 1.6];
        
        fs.n = fignum;
        fs.title = sprintf('%d - Effecs Increasing TF-VIIa',fignum);
%         fs.position = [926 245 447 329];
        fs.plots = {ps};
        fs.sub = [1, 1];
    end

main;
end

function ydot = odefun(~,y,p)
% Collect param values in cell array and redefine params with names
paramsCell = num2cell(p);
[k1,k2,k3,k4,k5,k6,k7,k8,k9,k10, ...
 k11,k12,k13,k14,k15,k16,k17,k18,k19,k20, ...
 km1,km1_2,km2,km3,km4,km5,km6,km7,km8]=paramsCell{:};

% dS;dPC;dAPC;dAPC_S;dAPC_S_V

% Collect y-vals in cell array and redefine y-vals with names
yCell = num2cell(y);
[TF_VIIa,IX,X,V,VIII,II,VIIIa_IXa,Va_Xa,IIa,Va_Xa_II,mIIa,...
 TF_VIIa_IX,TF_VIIa_X,VIIIa_IXa_X,IXa,Xa,Va,VIIIa,I,S,PC,APC,...
 APC_S,APC_S_V] = yCell{:};

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

% VIIIa-IXa - species 7 
if (k20) && (I < VIIIa_IXa)
    dVIIIa_IXa = 0;
else
    dVIIIa_IXa = k7*VIIIa*IXa - k9*VIIIa_IXa - k6*VIIIa_IXa*X ...
        + k18*VIIIa_IXa_X + k13*VIIIa_IXa_X ...
        + k20 * (-abs(I - VIIIa_IXa) + (I - VIIIa_IXa));
end

% Va-Xa (incorrect in original paper) - species 8
dVa_Xa = k8*Xa*Va - k10*Va_Xa + k19*Va_Xa_II ...
    - k6*Va_Xa*II + k14*Va_Xa_II;

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
dIXa = k9*VIIIa_IXa - k7*VIIIa*IXa + k11*TF_VIIa_IX ...
    + k15*(IX*Xa + IX*Va_Xa);

% Xa (incorrect in original paper: k8 should be used) - species 16
dXa = k10*Va_Xa - k8*Xa*Va + k12*TF_VIIa_X + k13*VIIIa_IXa_X;

% Va (incorrect in original paper: k8 should be used) - species 17
dVa = k10*Va_Xa - k8*Xa*Va + k1*V*Xa + k2*V*IIa + k2*V*mIIa;

% VIIIa - species 18
dVIIIa = k9*VIIIa_IXa - k7*VIIIa*IXa + k3*VIII*Xa ...
    + k4*(VIII*IIa + VIII*mIIa);

% maximal VIIIa_IXa
dI = k20 * (-abs(I - VIIIa_IXa) + (I - VIIIa_IXa));

%% MODIFICATIONS
% recycling upon APC complex dissosiation
dV = dV + km8 * APC_S_V;
dXa = dXa + km6 * Va_Xa * APC_S + km6 * APC_S_V * Va_Xa;

% Va Inactivation
dVa = dVa - km6 * APC_S * Va - km6 * APC_S_V * Va;
dVa_Xa = dVa_Xa - km6 * APC_S * Va_Xa - km6 * APC_S_V * Va_Xa;

% VIIa inactivation
dTF_VIIa = dTF_VIIa - km7 * APC_S_V * TF_VIIa;

% Protein C
dPC = -km1_2 * PC * mIIa - km1 * PC * IIa + km8 * (APC + APC_S + APC_S_V);

% active protein c
dAPC = km1_2 * PC * mIIa + km1 * PC * IIa - km2 * APC * S + km3 * APC_S ...
    - km8 * APC;

% active protein c - s
dAPC_S = km2 * APC * S - km4 * APC_S * V + km5 * APC_S_V - km8 * APC_S;

% active protein c - s - v
dAPC_S_V = km4 * APC_S * V - km5 * APC_S_V - km8 * APC_S_V;

% protein s
dS = -km2 * APC * S + km3 * APC_S + km8 * (APC_S + APC_S_V);

% Collect all ODEs for output
ydot = [dTF_VIIa;dIX;dX;dV;dVIII;dII;dVIIIa_IXa;dVa_Xa;dIIa;dVa_Xa_II; ...
        dmIIa;dTF_VIIa_IX;dTF_VIIa_X;dVIIIa_IXa_X;dIXa;dXa;dVa;dVIIIa; ...
        dI;dS;dPC;dAPC;dAPC_S;dAPC_S_V];
end

%% Figure Making Helper Functions
function fighand = makefigure(fs)
if isfield(fs, 'position')
    fighand = setupfig(fs.n, fs.title, fs.position);
else
    fighand = setupfig(fs.n, fs.title);
end
for i = 1:length(fs.plots)
    subplot(fs.sub(1), fs.sub(2), i);
    hold on;
    makeplot(fs.plots{i});
end
end

function makeplot(ps)
grid on;
if iscell(ps.x) && iscell(ps.y)
    n = length(ps.x);
    cycle = 1;
    spec = {'-'};
    fieldstr = 'spec_cycle';
    if isfield(ps, fieldstr)
        spec = ps.(fieldstr);
    end
    if  n ~= length(ps.y)
        error('x and y cell arrays have mismatched dimensions.');
    end
    for i = 1:n
        fieldstr = 'color_cycle';
        if isfield(ps, fieldstr)
            if (mod(i, ps.(fieldstr)) == 1) && (i ~= 1)
                ax = gca;
                ax.ColorOrderIndex = 1;
                cycle = cycle + 1;
            end
        end
            plot(ps.x{i}, ps.y{i},spec{cycle});
    end
else
    if iscell(ps.x) || iscell(ps.y)
        error('x and y values must both be matrices or cell arrays.');
    end
    plot(ps.x, ps.y);
end
xlabel(ps.xlabel);
ylabel(ps.ylabel);
xlim(ps.xlim);
ylim(ps.ylim);
leg = legend(ps.legend);
if isfield(ps, 'legend_title')
    title(leg, ps.legend_title);
end
leg.LineWidth = 0.5;
if isfield(ps, 'legend_loc')
    leg.Location = ps.legend_loc;
end
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
box off; grid on;
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

function corban_figure_defaults
% CORBANFIGUREDEFAULTS Sets default values to make pretty figures.
fontSize = 13;
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