%Provided function, not written by Shane Billingsley

% Interpolation GUI designed for Table A-6 (or any similar superheated
% table organized by pressure blocks, such as Table A-13.
%
% Run this file, (you defintely don't need to read it) and a window will
% appear that will interpolate properties for you. It is designed to look 
% like entries in Table A-6.
% 
% You can choose to interpolate between two different pressures, or for a 
% differnt property at a constant pressure. It will compute single 
% interpolations atfirst, but if you continue past that it will also do 
% double interpolations.
%
% Finally, this definitely works on the newer versions of MATLAB (I wrote
% this on 2023a), and I *think* it works on older versions, but I don't
% have the older version so I can't test it easily. If you have problems
% message me (Sarah) on Slack. To my knowledge, I'm not using any functions
% that were introduced after about 2011.

function Superheated_Table_Interpolation()
fig = uifigure("Name", 'Double Linear Interpolation', 'Position', [200 200 1100 500]);
mainGrid = uigridlayout(fig, 'ColumnWidth', {'1x', '1x', '1x', '1x', '1x', '1x'}, ...
    'RowHeight', {225, '1x'});

stepNumber = 1;

%% Properties Panel
prop_type_panel = uipanel(mainGrid, 'Title', 'Properties of Interest');
prop_type_panel.Layout.Row = 1;
prop_type_panel.Layout.Column = 3;
prop_type_grid = uigridlayout(prop_type_panel);
prop_type_grid.RowHeight = {22, 22, 22, 22};
prop_type_grid.ColumnWidth = {'1x'};

prop_nu = uicheckbox(prop_type_grid,'Text', 'Specific Volume', 'Value', 1, 'ValueChangedFcn', @(prop_nu, event) prop_changed());
prop_nu.Layout.Row = 1;
prop_u = uicheckbox(prop_type_grid,'Text', 'Internal Energy', 'Value', 1, 'ValueChangedFcn', @(prop_u, event) prop_changed());
prop_u.Layout.Row = 2;
prop_h = uicheckbox(prop_type_grid, 'Text', 'Enthalpy', 'Value', 1, 'ValueChangedFcn', @(prop_h, event) prop_changed());
prop_h.Layout.Row = 3;
prop_s = uicheckbox(prop_type_grid, 'Text', 'Entropy', 'Value', 1, 'ValueChangedFcn', @(prop_s, event) prop_changed());
prop_s.Layout.Row = 4;


%% First Type Panel
first_type_panel = uipanel(mainGrid, 'Title', 'First Interpolation Type');
first_type_panel.Layout.Row = 1;
first_type_panel.Layout.Column = 4;
first_type_grid = uigridlayout(first_type_panel);
first_type_grid.RowHeight = {84};
first_type_grid.ColumnWidth = {'1x'};

first_type_bg = uibuttongroup(first_type_grid, 'SelectionChangedFcn', @(first_type_bg, event) update());
first_type_bg.Layout.Row = 1;
first_type_bg.Layout.Column = 1;
first_type_bg.BorderType = "none";

first_type_rb_Pressure = uiradiobutton(first_type_bg,'Position',[10 60 120 15], 'Text', 'Cross-pressure');
first_type_rb_other = uiradiobutton(first_type_bg,'Position',[10 38 120 15], 'Text', 'Constant-pressure', 'Value', 1);


%% First Lookup Panel
first_lookup_panel = uipanel(mainGrid, 'Title', 'First Index Property');
first_lookup_panel.Layout.Row = 1;
first_lookup_panel.Layout.Column = 5;
first_lookup_grid = uigridlayout(first_lookup_panel);
first_lookup_grid.RowHeight = {162};
first_lookup_grid.ColumnWidth = {'1x'};

first_lookup_bg = uibuttongroup(first_lookup_grid, 'SelectionChangedFcn', @(first_lookup_bg, event) update());
first_lookup_bg.Layout.Row = 1;
first_lookup_bg.Layout.Column = 1;
first_lookup_bg.BorderType = "none";

first_lookup_rb_T = uiradiobutton(first_lookup_bg,'Position',[10 148 120 15], 'Text', 'Temperature');
first_lookup_rb_nu = uiradiobutton(first_lookup_bg,'Position',[10 126 120 15], 'Text', 'Specific Volume');
first_lookup_rb_u = uiradiobutton(first_lookup_bg,'Position',[10 104 120 15], 'Text', 'Internal Energy');
first_lookup_rb_h = uiradiobutton(first_lookup_bg,'Position',[10 82 120 15], 'Text', 'Enthalpy');
first_lookup_rb_s = uiradiobutton(first_lookup_bg,'Position',[10 60 120 15], 'Text', 'Entropy');
first_lookup_rb_P = uiradiobutton(first_lookup_bg,'Position',[10 38 120 15], 'Text', 'Pressure');


%% Second Lookup Panel
second_lookup_panel = uipanel(mainGrid, 'Title', 'Second Index Property');
second_lookup_panel.Layout.Row = 1;
second_lookup_panel.Layout.Column = 6;
second_lookup_grid = uigridlayout(second_lookup_panel);
second_lookup_grid.RowHeight = {162};
second_lookup_grid.ColumnWidth = {'1x'};

second_lookup_bg = uibuttongroup(second_lookup_grid, 'SelectionChangedFcn', @(second_lookup_bg, event) update());
second_lookup_bg.Layout.Row = 1;
second_lookup_bg.Layout.Column = 1;
second_lookup_bg.BorderType = "none";

second_lookup_rb_T = uiradiobutton(second_lookup_bg,'Position',[10 148 120 15], 'Text', 'Temperature');
second_lookup_rb_nu = uiradiobutton(second_lookup_bg,'Position',[10 126 120 15], 'Text', 'Specific Volume');
second_lookup_rb_u = uiradiobutton(second_lookup_bg,'Position',[10 104 120 15], 'Text', 'Internal Energy');
second_lookup_rb_h = uiradiobutton(second_lookup_bg,'Position',[10 82 120 15], 'Text', 'Enthalpy');
second_lookup_rb_s = uiradiobutton(second_lookup_bg,'Position',[10 60 120 15], 'Text', 'Entropy');
second_lookup_rb_P = uiradiobutton(second_lookup_bg,'Position',[10 38 120 15], 'Text', 'Pressure');


%% Instruction Detail Panel

detail_panel = uipanel(mainGrid, 'Title', 'Instructions');
detail_panel.Layout.Row = 1;
detail_panel.Layout.Column = [1 2];
detail_grid = uigridlayout(detail_panel, 'RowSpacing', 0);
detail_grid.RowHeight = {22, 2, 22, 22, 22, 22, 22, '1x', 30};
detail_grid.ColumnWidth = {'1x', '1x'};

detail_description_label_Title = uilabel(detail_grid, 'FontWeight', 'bold');
detail_description_label_Title.Layout.Row = 1;
detail_description_label_Title.Layout.Column = [1 2];

detail_description_label_1 = uilabel(detail_grid);
detail_description_label_1.Layout.Row = 3;
detail_description_label_1.Layout.Column = [1 2];

detail_description_label_2 = uilabel(detail_grid);
detail_description_label_2.Layout.Row = 4;
detail_description_label_2.Layout.Column = [1 2];

detail_description_label_3 = uilabel(detail_grid);
detail_description_label_3.Layout.Row = 5;
detail_description_label_3.Layout.Column = [1 2];

detail_description_label_4 = uilabel(detail_grid);
detail_description_label_4.Layout.Row = 6;
detail_description_label_4.Layout.Column = [1 2];

detail_description_label_5 = uilabel(detail_grid);
detail_description_label_5.Layout.Row = 7;
detail_description_label_5.Layout.Column = [1 2];

prev_button = uibutton(detail_grid, 'push', 'Text', 'Go back', 'ButtonPushedFcn', @(event, prev_button) back(), 'Visible', 'off');
prev_button.Layout.Row = 9;
prev_button.Layout.Column = 1;

next_button = uibutton(detail_grid, 'push', 'Text', 'Next', 'ButtonPushedFcn', @(event, next_button) advance());
next_button.Layout.Row = 9;
next_button.Layout.Column = 2;

stepChange()
%% P1 Panel
P1_panel = uipanel(mainGrid);
P1_panel.Layout.Row = 2;
P1_panel.Layout.Column = [1 2];

P1_panel_grid = uigridlayout(P1_panel);
P1_panel_grid.RowHeight = {20, 40, 20, 20, 20, 20, '1x'};
P1_panel_grid.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};

P1_label_grid = uigridlayout(P1_panel_grid);
P1_label_grid.Layout.Row = 1;
P1_label_grid.Layout.Column = [1 5];
P1_label_grid.RowHeight = {'1x'};
P1_label_grid.ColumnWidth = {'1x', 50, '1x'};
P1_label_grid.Padding = [0 0 0 0];

P1_label = uilabel(P1_label_grid, 'Text', 'P =', 'HorizontalAlignment', 'right');
P1_label.Layout.Column = 1;

P1_pressure = uieditfield(P1_label_grid, 'numeric', 'Value', 1.00, "ValueDisplayFormat","%g", 'ValueChangedFcn', @(event, P1_pressure) update());
P1_pressure.Layout.Column = 2;

P1_unit_label = uilabel(P1_label_grid, 'Text', 'MPa');
P1_unit_label.Layout.Column = 3;


P1_T_label = uilabel(P1_panel_grid, 'Text', {'T'; append(char(176), 'C')}, ...
    'HorizontalAlignment', 'left');
P1_T_label.Layout.Row = 2;
P1_T_label.Layout.Column = 1;

P1_nu_label = uilabel(P1_panel_grid, 'Text', {char(957); append('m', char(179), '/kg')}, ...
    'HorizontalAlignment', 'left');
P1_nu_label.Layout.Row = 2;
P1_nu_label.Layout.Column = 2;

P1_u_label = uilabel(P1_panel_grid, 'Text', {'u'; 'kJ/kg'}, ...
    'HorizontalAlignment', 'left');
P1_u_label.Layout.Row = 2;
P1_u_label.Layout.Column = 3;

P1_h_label = uilabel(P1_panel_grid, 'Text', {'h'; 'kJ/kg'}, ...
    'HorizontalAlignment', 'left');
P1_h_label.Layout.Row = 2;
P1_h_label.Layout.Column = 4;

P1_s_label = uilabel(P1_panel_grid, 'Text', {'s'; append('kJ/kg', char(183), 'K')}, ...
    'HorizontalAlignment', 'left');
P1_s_label.Layout.Row = 2;
P1_s_label.Layout.Column = 5;


P1_T1_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 200,  'ValueChangedFcn', @(event, P1_T1_field) update(), "ValueDisplayFormat","%g");
P1_T1_field.Layout.Row = 3;
P1_T1_field.Layout.Column = 1;

P1_T2_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 250,  'ValueChangedFcn', @(event, P1_T2_field) update(), "ValueDisplayFormat","%g");
P1_T2_field.Layout.Row = 5;
P1_T2_field.Layout.Column = 1;

P1_Tx_field = uieditfield(P1_panel_grid, 'numeric', 'Value', (P1_T1_field.Value + P1_T2_field.Value)/2,  'ValueChangedFcn', @(event, P1_Tx_field) update(), "ValueDisplayFormat","%g");
P1_Tx_field.Layout.Row = 4;
P1_Tx_field.Layout.Column = 1;

P1_nu1_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 0.20602,  'ValueChangedFcn', @(event, P1_nu1_field) update(), "ValueDisplayFormat","%.5f");
P1_nu1_field.Layout.Row = 3;
P1_nu1_field.Layout.Column = 2;

P1_nu2_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 0.23275,  'ValueChangedFcn', @(event, P1_nu2_field) update(), "ValueDisplayFormat","%.5f");
P1_nu2_field.Layout.Row = 5;
P1_nu2_field.Layout.Column = 2;

P1_nux_field = uieditfield(P1_panel_grid, 'numeric', 'Value', (P1_nu1_field.Value + P1_nu2_field.Value)/2,  'ValueChangedFcn', @(event, P1_nux_field) update(), "ValueDisplayFormat","%.5f");
P1_nux_field.Layout.Row = 4;
P1_nux_field.Layout.Column = 2;

P1_u1_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 2622.3,  'ValueChangedFcn', @(event, P1_u1_field) update(), "ValueDisplayFormat","%.1f");
P1_u1_field.Layout.Row = 3;
P1_u1_field.Layout.Column = 3;

P1_u2_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 2710.4,  'ValueChangedFcn', @(event, P1_u2_field) update(), "ValueDisplayFormat","%.1f");
P1_u2_field.Layout.Row = 5;
P1_u2_field.Layout.Column = 3;

P1_ux_field = uieditfield(P1_panel_grid, 'numeric', 'Value', (P1_u1_field.Value + P1_u2_field.Value)/2,  'ValueChangedFcn', @(event, P1_ux_field) update(), "ValueDisplayFormat","%.1f");
P1_ux_field.Layout.Row = 4;
P1_ux_field.Layout.Column = 3;

P1_h1_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 2828.3,  'ValueChangedFcn', @(event, P1_h1_field) update(), "ValueDisplayFormat","%.1f");
P1_h1_field.Layout.Row = 3;
P1_h1_field.Layout.Column = 4;

P1_h2_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 2943.1,  'ValueChangedFcn', @(event, P1_h2_field) update(), "ValueDisplayFormat","%.1f");
P1_h2_field.Layout.Row = 5;
P1_h2_field.Layout.Column = 4;

P1_hx_field = uieditfield(P1_panel_grid, 'numeric', 'Value', (P1_h1_field.Value + P1_h2_field.Value)/2,  'ValueChangedFcn', @(event, P1_hx_field) update(), "ValueDisplayFormat","%.1f");
P1_hx_field.Layout.Row = 4;
P1_hx_field.Layout.Column = 4;

P1_s1_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 6.6956,  'ValueChangedFcn', @(event, P1_s1_field) update(), "ValueDisplayFormat","%.4f");
P1_s1_field.Layout.Row = 3;
P1_s1_field.Layout.Column = 5;

P1_s2_field = uieditfield(P1_panel_grid, 'numeric', 'Value', 6.9265,  'ValueChangedFcn', @(event, P1_s2_field) update(), "ValueDisplayFormat","%.4f");
P1_s2_field.Layout.Row = 5;
P1_s2_field.Layout.Column = 5;

P1_sx_field = uieditfield(P1_panel_grid, 'numeric', 'Value', (P1_s1_field.Value + P1_s2_field.Value)/2,  'ValueChangedFcn', @(event, P1_sx_field) update(), "ValueDisplayFormat","%.4f");
P1_sx_field.Layout.Row = 4;
P1_sx_field.Layout.Column = 5;


%% P2 Panel
P2_panel = uipanel(mainGrid);
P2_panel.Layout.Row = 2;
P2_panel.Layout.Column = [5 6];

P2_panel_grid = uigridlayout(P2_panel);
P2_panel_grid.RowHeight = {20, 40, 20, 20, 20, '1x'};
P2_panel_grid.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};

P2_label_grid = uigridlayout(P2_panel_grid);
P2_label_grid.Layout.Row = 1;
P2_label_grid.Layout.Column = [1 5];
P2_label_grid.RowHeight = {'1x'};
P2_label_grid.ColumnWidth = {'1x', 50, '1x'};
P2_label_grid.Padding = [0 0 0 0];

P2_label = uilabel(P2_label_grid, 'Text', 'P =', 'HorizontalAlignment', 'right');
P2_label.Layout.Column = 1;

P2_pressure = uieditfield(P2_label_grid, 'numeric', 'Value', 1.20, "ValueDisplayFormat","%g",  'ValueChangedFcn', @(event, P2_pressure) update());
P2_pressure.Layout.Column = 2;

P2_unit_label = uilabel(P2_label_grid, 'Text', 'MPa');
P2_unit_label.Layout.Column = 3;


P2_T_label = uilabel(P2_panel_grid, 'Text', {'T'; append(char(176), 'C')}, ...
    'HorizontalAlignment', 'left');
P2_T_label.Layout.Row = 2;
P2_T_label.Layout.Column = 1;

P2_nu_label = uilabel(P2_panel_grid, 'Text', {char(957); append('m', char(179), '/kg')}, ...
    'HorizontalAlignment', 'left');
P2_nu_label.Layout.Row = 2;
P2_nu_label.Layout.Column = 2;

P2_u_label = uilabel(P2_panel_grid, 'Text', {'u'; 'kJ/kg'}, ...
    'HorizontalAlignment', 'left');
P2_u_label.Layout.Row = 2;
P2_u_label.Layout.Column = 3;

P2_h_label = uilabel(P2_panel_grid, 'Text', {'h'; 'kJ/kg'}, ...
    'HorizontalAlignment', 'left');
P2_h_label.Layout.Row = 2;
P2_h_label.Layout.Column = 4;

P2_s_label = uilabel(P2_panel_grid, 'Text', {'s'; append('kJ/kg', char(183), 'K')}, ...
    'HorizontalAlignment', 'left');
P2_s_label.Layout.Row = 2;
P2_s_label.Layout.Column = 5;


P2_T1_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 200,  'ValueChangedFcn', @(event, P2_T1_field) update(), "ValueDisplayFormat","%g");
P2_T1_field.Layout.Row = 3;
P2_T1_field.Layout.Column = 1;

P2_T2_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 250,  'ValueChangedFcn', @(event, P2_T2_field) update(), "ValueDisplayFormat","%g");
P2_T2_field.Layout.Row = 5;
P2_T2_field.Layout.Column = 1;

P2_Tx_field = uieditfield(P2_panel_grid, 'numeric', 'Value', (P2_T1_field.Value + P2_T2_field.Value)/2,  'ValueChangedFcn', @(event, P2_Tx_field) update(), "ValueDisplayFormat","%g");
P2_Tx_field.Layout.Row = 4;
P2_Tx_field.Layout.Column = 1;

P2_nu1_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 0.16934,  'ValueChangedFcn', @(event, P2_nu1_field) update(), "ValueDisplayFormat","%.5f");
P2_nu1_field.Layout.Row = 3;
P2_nu1_field.Layout.Column = 2;

P2_nu2_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 0.19241,  'ValueChangedFcn', @(event, P2_nu2_field) update(), "ValueDisplayFormat","%.5f");
P2_nu2_field.Layout.Row = 5;
P2_nu2_field.Layout.Column = 2;

P2_nux_field = uieditfield(P2_panel_grid, 'numeric', 'Value', (P2_nu1_field.Value + P2_nu2_field.Value)/2,  'ValueChangedFcn', @(event, P2_nux_field) update(), "ValueDisplayFormat","%.5f");
P2_nux_field.Layout.Row = 4;
P2_nux_field.Layout.Column = 2;

P2_u1_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 2612.9,  'ValueChangedFcn', @(event, P2_u1_field) update(), "ValueDisplayFormat","%.1f");
P2_u1_field.Layout.Row = 3;
P2_u1_field.Layout.Column = 3;

P2_u2_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 2704.7,  'ValueChangedFcn', @(event, P2_u2_field) update(), "ValueDisplayFormat","%.1f");
P2_u2_field.Layout.Row = 5;
P2_u2_field.Layout.Column = 3;

P2_ux_field = uieditfield(P2_panel_grid, 'numeric', 'Value', (P2_u1_field.Value + P2_u2_field.Value)/2,  'ValueChangedFcn', @(event, P2_ux_field) update(), "ValueDisplayFormat","%.1f");
P2_ux_field.Layout.Row = 4;
P2_ux_field.Layout.Column = 3;

P2_h1_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 2816.1,  'ValueChangedFcn', @(event, P2_h1_field) update(), "ValueDisplayFormat","%.1f");
P2_h1_field.Layout.Row = 3;
P2_h1_field.Layout.Column = 4;

P2_h2_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 2935.6,  'ValueChangedFcn', @(event, P2_h1_field) update(), "ValueDisplayFormat","%.1f");
P2_h2_field.Layout.Row = 5;
P2_h2_field.Layout.Column = 4;

P2_hx_field = uieditfield(P2_panel_grid, 'numeric', 'Value', (P2_h1_field.Value + P2_h2_field.Value)/2,  'ValueChangedFcn', @(event, P2_h1_field) update(), "ValueDisplayFormat","%.1f");
P2_hx_field.Layout.Row = 4;
P2_hx_field.Layout.Column = 4;

P2_s1_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 6.5909,  'ValueChangedFcn', @(event, P2_s1_field) update(), "ValueDisplayFormat","%.4f");
P2_s1_field.Layout.Row = 3;
P2_s1_field.Layout.Column = 5;

P2_s2_field = uieditfield(P2_panel_grid, 'numeric', 'Value', 6.8313,  'ValueChangedFcn', @(event, P2_s2_field) update(), "ValueDisplayFormat","%.4f");
P2_s2_field.Layout.Row = 5;
P2_s2_field.Layout.Column = 5;

P2_sx_field = uieditfield(P2_panel_grid, 'numeric', 'Value', (P2_s1_field.Value + P2_s2_field.Value)/2,  'ValueChangedFcn', @(event, P2_sx_field) update(), "ValueDisplayFormat","%.4f");
P2_sx_field.Layout.Row = 4;
P2_sx_field.Layout.Column = 5;


%% Px Panel
Px_panel = uipanel(mainGrid);
Px_panel.Layout.Row = 2;
Px_panel.Layout.Column = [3 4];

Px_panel_grid = uigridlayout(Px_panel);
Px_panel_grid.RowHeight = {20, 40, 20, 20, 20, '1x'};
Px_panel_grid.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};

Px_label_grid = uigridlayout(Px_panel_grid);
Px_label_grid.Layout.Row = 1;
Px_label_grid.Layout.Column = [1 5];
Px_label_grid.RowHeight = {'1x'};
Px_label_grid.ColumnWidth = {'1x', 50, '1x'};
Px_label_grid.Padding = [0 0 0 0];

Px_label = uilabel(Px_label_grid, 'Text', 'P =', 'HorizontalAlignment', 'right');
Px_label.Layout.Column = 1;

Px_pressure = uieditfield(Px_label_grid, 'numeric', 'Value', 1.20, "ValueDisplayFormat","%g",  'ValueChangedFcn', @(event, Px_pressure) update());
Px_pressure.Layout.Column = 2;

Px_unit_label = uilabel(Px_label_grid, 'Text', 'MPa');
Px_unit_label.Layout.Column = 3;


Px_T_label = uilabel(Px_panel_grid, 'Text', {'T'; append(char(176), 'C')}, ...
    'HorizontalAlignment', 'left');
Px_T_label.Layout.Row = 2;
Px_T_label.Layout.Column = 1;

Px_nu_label = uilabel(Px_panel_grid, 'Text', {char(957); append('m', char(179), '/kg')}, ...
    'HorizontalAlignment', 'left');
Px_nu_label.Layout.Row = 2;
Px_nu_label.Layout.Column = 2;

Px_u_label = uilabel(Px_panel_grid, 'Text', {'u'; 'kJ/kg'}, ...
    'HorizontalAlignment', 'left');
Px_u_label.Layout.Row = 2;
Px_u_label.Layout.Column = 3;

Px_h_label = uilabel(Px_panel_grid, 'Text', {'h'; 'kJ/kg'}, ...
    'HorizontalAlignment', 'left');
Px_h_label.Layout.Row = 2;
Px_h_label.Layout.Column = 4;

Px_s_label = uilabel(Px_panel_grid, 'Text', {'s'; append('kJ/kg', char(183), 'K')}, ...
    'HorizontalAlignment', 'left');
Px_s_label.Layout.Row = 2;
Px_s_label.Layout.Column = 5;


Px_T1_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 200,  'ValueChangedFcn', @(event, Px_T1_field) update(), "ValueDisplayFormat","%g");
Px_T1_field.Layout.Row = 3;
Px_T1_field.Layout.Column = 1;

Px_T2_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 250,  'ValueChangedFcn', @(event, Px_T2_field) update(), "ValueDisplayFormat","%g");
Px_T2_field.Layout.Row = 5;
Px_T2_field.Layout.Column = 1;

Px_Tx_field = uieditfield(Px_panel_grid, 'numeric', 'Value', (Px_T1_field.Value + Px_T2_field.Value)/2,  'ValueChangedFcn', @(event, Px_Tx_field) update(), "ValueDisplayFormat","%g");
Px_Tx_field.Layout.Row = 4;
Px_Tx_field.Layout.Column = 1;

Px_nu1_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 0.16934,  'ValueChangedFcn', @(event, Px_nu1_field) update(), "ValueDisplayFormat","%.5f");
Px_nu1_field.Layout.Row = 3;
Px_nu1_field.Layout.Column = 2;

Px_nu2_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 0.19241,  'ValueChangedFcn', @(event, Px_nu2_field) update(), "ValueDisplayFormat","%.5f");
Px_nu2_field.Layout.Row = 5;
Px_nu2_field.Layout.Column = 2;

Px_nux_field = uieditfield(Px_panel_grid, 'numeric', 'Value', (Px_nu1_field.Value + Px_nu2_field.Value)/2,  'ValueChangedFcn', @(event, Px_nux_field) update(), "ValueDisplayFormat","%.5f");
Px_nux_field.Layout.Row = 4;
Px_nux_field.Layout.Column = 2;

Px_u1_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 2612.9,  'ValueChangedFcn', @(event, Px_u1_field) update(), "ValueDisplayFormat","%.1f");
Px_u1_field.Layout.Row = 3;
Px_u1_field.Layout.Column = 3;

Px_u2_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 2704.7,  'ValueChangedFcn', @(event, Px_u2_field) update(), "ValueDisplayFormat","%.1f");
Px_u2_field.Layout.Row = 5;
Px_u2_field.Layout.Column = 3;

Px_ux_field = uieditfield(Px_panel_grid, 'numeric', 'Value', (Px_u1_field.Value + Px_u2_field.Value)/2,  'ValueChangedFcn', @(event, Px_ux_field) update(), "ValueDisplayFormat","%.1f");
Px_ux_field.Layout.Row = 4;
Px_ux_field.Layout.Column = 3;

Px_h1_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 2816.1,  'ValueChangedFcn', @(event, Px_h1_field) update(), "ValueDisplayFormat","%.1f");
Px_h1_field.Layout.Row = 3;
Px_h1_field.Layout.Column = 4;

Px_h2_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 2935.6,  'ValueChangedFcn', @(event, Px_h2_field) update(), "ValueDisplayFormat","%.1f");
Px_h2_field.Layout.Row = 5;
Px_h2_field.Layout.Column = 4;

Px_hx_field = uieditfield(Px_panel_grid, 'numeric', 'Value', (Px_h1_field.Value + Px_h2_field.Value)/2,  'ValueChangedFcn', @(event, Px_hx_field) update(), "ValueDisplayFormat","%.1f");
Px_hx_field.Layout.Row = 4;
Px_hx_field.Layout.Column = 4;

Px_s1_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 6.5909,  'ValueChangedFcn', @(event, Px_s1_field) update(), "ValueDisplayFormat","%.4f");
Px_s1_field.Layout.Row = 3;
Px_s1_field.Layout.Column = 5;

Px_s2_field = uieditfield(Px_panel_grid, 'numeric', 'Value', 6.8313,  'ValueChangedFcn', @(event, Px_s2_field) update(), "ValueDisplayFormat","%.4f");
Px_s2_field.Layout.Row = 5;
Px_s2_field.Layout.Column = 5;

Px_sx_field = uieditfield(Px_panel_grid, 'numeric', 'Value', (Px_s1_field.Value + Px_s2_field.Value)/2,  'ValueChangedFcn', @(event, Px_sx_field) update(), "ValueDisplayFormat","%.4f");
Px_sx_field.Layout.Row = 4;
Px_sx_field.Layout.Column = 5;

update()

    function prop_changed()
        stepNumber = 1;
        update();
    end

    function stepChange()
        next_button.Enable = 'off';
        prev_button.Enable = 'off';
        if stepNumber == 1
            detail_description_label_Title.Text = 'Step 1:';
            
            detail_description_label_1.Text = 'Select the properties you are interested in finding, including';
            detail_description_label_2.Text = 'the one(s) you need to index by. You can leave these all';
            detail_description_label_3.Text = 'selected if you are willing to enter all of them.';
            detail_description_label_4.Text = '';
            detail_description_label_5.Text = '';
        elseif stepNumber == 2
            detail_description_label_Title.Text = 'Step 2:';
            
            detail_description_label_1.Text = 'Choose the layout of the intrpolation you would like to';
            detail_description_label_2.Text = 'compute. Cross-pressure (between different pressure blocks),';
            detail_description_label_3.Text = 'or constant-pressure. Enter the surrounding properties.';
            detail_description_label_4.Text = '';
            detail_description_label_5.Text = '';
        elseif stepNumber == 3
            detail_description_label_Title.Text = 'Step 3:';
            
            detail_description_label_1.Text = 'Choose the property you want to index by. Enter this property,';
            detail_description_label_2.Text = 'and the first interpolation values will populate as soon as you';
            detail_description_label_3.Text = 'click outside the box.';
            detail_description_label_4.Text = '';
            detail_description_label_5.Text = 'If you want to continue to a double-interpolation, click Next.';
        else
            detail_description_label_Title.Text = 'Step 4:';

            detail_description_label_1.Text = 'Select the second index property.';
            detail_description_label_2.Text = '';
            detail_description_label_3.Text = 'Enter the remaining surrounding properties, and the value of';
            detail_description_label_4.Text = 'the second index property. The interpolated values will';
            detail_description_label_5.Text = 'calculate automatically.';
        end
        update();
        pause(0.1);
        next_button.Enable = 'on';
        prev_button.Enable = 'on';
    end

    function advance()
        stepNumber = stepNumber + 1;
        if stepNumber == 4
            next_button.Visible = 'off';
            prev_button.Visible = 'on';
        elseif stepNumber == 1
            next_button.Visible = 'on';
            prev_button.Visible = 'off';
        else
            next_button.Visible = 'on';
            prev_button.Visible = 'on';
        end
        stepChange()
    end

    function back()
        stepNumber = stepNumber - 1;
        if stepNumber == 6
            next_button.Visible = 'off';
            prev_button.Visible = 'on';
        elseif stepNumber == 1
            next_button.Visible = 'on';
            prev_button.Visible = 'off';
        else
            next_button.Visible = 'on';
            prev_button.Visible = 'on';
        end
        stepChange()
    end

    function yq = interpolate(x1, x2, y1, y2, xq)
        yq = y1 + (xq - x1)*(y2 - y1)/(x2 - x1);
    end

    function update()
        if stepNumber == 1
            P2_panel.Visible = 'off';
            Px_panel.Visible = 'off';

            prop_type_panel.Visible = 'on';
            prop_type_panel.Enable = 'on';
            first_type_panel.Visible = 'off';
            first_lookup_panel.Visible = 'off';
            second_lookup_panel.Visible = 'off';

            P1_Tx_field.Enable = 'off';
            P1_nux_field.Enable = 'off';
            P1_ux_field.Enable = 'off';
            P1_hx_field.Enable = 'off';
            P1_sx_field.Enable = 'off';

            P2_Tx_field.Enable = 'off';
            P2_nux_field.Enable = 'off';
            P2_ux_field.Enable = 'off';
            P2_hx_field.Enable = 'off';
            P2_sx_field.Enable = 'off';

            Px_Tx_field.Enable = 'off';
            Px_nux_field.Enable = 'off';
            Px_ux_field.Enable = 'off';
            Px_hx_field.Enable = 'off';
            Px_sx_field.Enable = 'off';

            if (prop_nu.Value + prop_u.Value + prop_h.Value + prop_s.Value) == 0
                prop_h.Value = 1;
            end

            if prop_nu.Value == 0
                P1_nu1_field.Visible = 'off';
                P1_nux_field.Visible = 'off';
                P1_nu2_field.Visible = 'off';
                P2_nu1_field.Visible = 'off';
                P2_nux_field.Visible = 'off';
                P2_nu2_field.Visible = 'off';
                Px_nu1_field.Visible = 'off';
                Px_nux_field.Visible = 'off';
                Px_nu2_field.Visible = 'off';
            else
                P1_nu1_field.Visible = 'on';
                P1_nux_field.Visible = 'on';
                P1_nu2_field.Visible = 'on';
                P2_nu1_field.Visible = 'on';
                P2_nux_field.Visible = 'on';
                P2_nu2_field.Visible = 'on';
                Px_nu1_field.Visible = 'on';
                Px_nux_field.Visible = 'on';
                Px_nu2_field.Visible = 'on';
            end

            if prop_u.Value == 0
                P1_u1_field.Visible = 'off';
                P1_ux_field.Visible = 'off';
                P1_u2_field.Visible = 'off';
                P2_u1_field.Visible = 'off';
                P2_ux_field.Visible = 'off';
                P2_u2_field.Visible = 'off';
                Px_u1_field.Visible = 'off';
                Px_ux_field.Visible = 'off';
                Px_u2_field.Visible = 'off';
            else
                P1_u1_field.Visible = 'on';
                P1_ux_field.Visible = 'on';
                P1_u2_field.Visible = 'on';
                P2_u1_field.Visible = 'on';
                P2_ux_field.Visible = 'on';
                P2_u2_field.Visible = 'on';
                Px_u1_field.Visible = 'on';
                Px_ux_field.Visible = 'on';
                Px_u2_field.Visible = 'on';
            end

            if prop_h.Value == 0
                P1_h1_field.Visible = 'off';
                P1_hx_field.Visible = 'off';
                P1_h2_field.Visible = 'off';
                P2_h1_field.Visible = 'off';
                P2_hx_field.Visible = 'off';
                P2_h2_field.Visible = 'off';
                Px_h1_field.Visible = 'off';
                Px_hx_field.Visible = 'off';
                Px_h2_field.Visible = 'off';
            else
                P1_h1_field.Visible = 'on';
                P1_hx_field.Visible = 'on';
                P1_h2_field.Visible = 'on';
                P2_h1_field.Visible = 'on';
                P2_hx_field.Visible = 'on';
                P2_h2_field.Visible = 'on';
                Px_h1_field.Visible = 'on';
                Px_hx_field.Visible = 'on';
                Px_h2_field.Visible = 'on';
            end

            if prop_s.Value == 0
                P1_s1_field.Visible = 'off';
                P1_sx_field.Visible = 'off';
                P1_s2_field.Visible = 'off';
                P2_s1_field.Visible = 'off';
                P2_sx_field.Visible = 'off';
                P2_s2_field.Visible = 'off';
                Px_s1_field.Visible = 'off';
                Px_sx_field.Visible = 'off';
                Px_s2_field.Visible = 'off';
            else
                P1_s1_field.Visible = 'on';
                P1_sx_field.Visible = 'on';
                P1_s2_field.Visible = 'on';
                P2_s1_field.Visible = 'on';
                P2_sx_field.Visible = 'on';
                P2_s2_field.Visible = 'on';
                Px_s1_field.Visible = 'on';
                Px_sx_field.Visible = 'on';
                Px_s2_field.Visible = 'on';
            end
            P1_Tx_field.Visible = 'on';
            P1_T2_field.Visible = 'on';
            P2_Tx_field.Visible = 'on';
            P2_T2_field.Visible = 'on';
            Px_Tx_field.Visible = 'on';
            Px_T2_field.Visible = 'on';
        elseif stepNumber == 2
            prop_type_panel.Visible = 'on';
            prop_type_panel.Enable = 'off';
            first_type_panel.Visible = 'on';
            first_type_panel.Enable = 'on';
            first_lookup_panel.Visible = 'off';
            second_lookup_panel.Visible = 'off';

            if first_type_rb_other.Value
                first_lookup_rb_P.Visible = 'off';
                first_lookup_rb_T.Visible = 'on';
                if first_lookup_rb_P.Value
                    first_lookup_rb_T.Value = 1;
                end
                P2_panel.Visible = 'off';
                Px_panel.Visible = 'off';

                P1_Tx_field.Enable = 'off';
                P1_nux_field.Enable = 'off';
                P1_ux_field.Enable = 'off';
                P1_hx_field.Enable = 'off';
                P1_sx_field.Enable = 'off';

                P2_Tx_field.Enable = 'off';
                P2_nux_field.Enable = 'off';
                P2_ux_field.Enable = 'off';
                P2_hx_field.Enable = 'off';
                P2_sx_field.Enable = 'off';

                Px_Tx_field.Enable = 'off';
                Px_nux_field.Enable = 'off';
                Px_ux_field.Enable = 'off';
                Px_hx_field.Enable = 'off';
                Px_sx_field.Enable = 'off';

                if (prop_nu.Value + prop_u.Value + prop_h.Value + prop_s.Value) == 0
                    prop_h.Value = 1;
                end

                first_lookup_rb_h.Visible = 'on';
                first_lookup_rb_nu.Visible = 'on';
                first_lookup_rb_u.Visible = 'on';
                first_lookup_rb_h.Visible = 'on';
                first_lookup_rb_s.Visible = 'on';

                if prop_nu.Value == 0

                    P1_nu1_field.Visible = 'off';
                    P1_nux_field.Visible = 'off';
                    P1_nu2_field.Visible = 'off';
                    P2_nu1_field.Visible = 'off';
                    P2_nux_field.Visible = 'off';
                    P2_nu2_field.Visible = 'off';
                    Px_nu1_field.Visible = 'off';
                    Px_nux_field.Visible = 'off';
                    Px_nu2_field.Visible = 'off';
                else
                    P1_nu1_field.Visible = 'on';
                    P1_nux_field.Visible = 'on';
                    P1_nu2_field.Visible = 'on';
                    P2_nu1_field.Visible = 'on';
                    P2_nux_field.Visible = 'on';
                    P2_nu2_field.Visible = 'on';
                    Px_nu1_field.Visible = 'on';
                    Px_nux_field.Visible = 'on';
                    Px_nu2_field.Visible = 'on';
                end

                if prop_u.Value == 0

                    P1_u1_field.Visible = 'off';
                    P1_ux_field.Visible = 'off';
                    P1_u2_field.Visible = 'off';
                    P2_u1_field.Visible = 'off';
                    P2_ux_field.Visible = 'off';
                    P2_u2_field.Visible = 'off';
                    Px_u1_field.Visible = 'off';
                    Px_ux_field.Visible = 'off';
                    Px_u2_field.Visible = 'off';
                else
                    P1_u1_field.Visible = 'on';
                    P1_ux_field.Visible = 'on';
                    P1_u2_field.Visible = 'on';
                    P2_u1_field.Visible = 'on';
                    P2_ux_field.Visible = 'on';
                    P2_u2_field.Visible = 'on';
                    Px_u1_field.Visible = 'on';
                    Px_ux_field.Visible = 'on';
                    Px_u2_field.Visible = 'on';
                end

                if prop_h.Value == 0

                    P1_h1_field.Visible = 'off';
                    P1_hx_field.Visible = 'off';
                    P1_h2_field.Visible = 'off';
                    P2_h1_field.Visible = 'off';
                    P2_hx_field.Visible = 'off';
                    P2_h2_field.Visible = 'off';
                    Px_h1_field.Visible = 'off';
                    Px_hx_field.Visible = 'off';
                    Px_h2_field.Visible = 'off';
                else
                    P1_h1_field.Visible = 'on';
                    P1_hx_field.Visible = 'on';
                    P1_h2_field.Visible = 'on';
                    P2_h1_field.Visible = 'on';
                    P2_hx_field.Visible = 'on';
                    P2_h2_field.Visible = 'on';
                    Px_h1_field.Visible = 'on';
                    Px_hx_field.Visible = 'on';
                    Px_h2_field.Visible = 'on';
                end

                if prop_s.Value == 0

                    P1_s1_field.Visible = 'off';
                    P1_sx_field.Visible = 'off';
                    P1_s2_field.Visible = 'off';
                    P2_s1_field.Visible = 'off';
                    P2_sx_field.Visible = 'off';
                    P2_s2_field.Visible = 'off';
                    Px_s1_field.Visible = 'off';
                    Px_sx_field.Visible = 'off';
                    Px_s2_field.Visible = 'off';
                else
                    P1_s1_field.Visible = 'on';
                    P1_sx_field.Visible = 'on';
                    P1_s2_field.Visible = 'on';
                    P2_s1_field.Visible = 'on';
                    P2_sx_field.Visible = 'on';
                    P2_s2_field.Visible = 'on';
                    Px_s1_field.Visible = 'on';
                    Px_sx_field.Visible = 'on';
                    Px_s2_field.Visible = 'on';
                end
                P1_Tx_field.Visible = 'on';
                P1_T2_field.Visible = 'on';
                P2_Tx_field.Visible = 'on';
                P2_T2_field.Visible = 'on';
                Px_Tx_field.Visible = 'on';
                Px_T2_field.Visible = 'on';

                Px_T1_field.Enable = 'on';
                Px_nu1_field.Enable = 'on';
                Px_u1_field.Enable = 'on';
                Px_h1_field.Enable = 'on';
                Px_s1_field.Enable = 'on';
            else
                first_lookup_rb_T.Visible = 'off';
                first_lookup_rb_P.Visible = 'on';
                if first_lookup_rb_T.Value
                    first_lookup_rb_P.Value = 1;
                end

                P2_panel.Visible = 'on';
                Px_panel.Visible = 'on';

                P1_Tx_field.Visible = 'off';
                P1_T2_field.Visible = 'off';
                P2_Tx_field.Visible = 'off';
                P2_T2_field.Visible = 'off';
                Px_Tx_field.Visible = 'off';
                Px_T2_field.Visible = 'off';

                P1_nux_field.Visible = 'off';
                P1_nu2_field.Visible = 'off';
                P2_nux_field.Visible = 'off';
                P2_nu2_field.Visible = 'off';
                Px_nux_field.Visible = 'off';
                Px_nu2_field.Visible = 'off';

                P1_ux_field.Visible = 'off';
                P1_u2_field.Visible = 'off';
                P2_ux_field.Visible = 'off';
                P2_u2_field.Visible = 'off';
                Px_ux_field.Visible = 'off';
                Px_u2_field.Visible = 'off';

                P1_hx_field.Visible = 'off';
                P1_h2_field.Visible = 'off';
                P2_hx_field.Visible = 'off';
                P2_h2_field.Visible = 'off';
                Px_hx_field.Visible = 'off';
                Px_h2_field.Visible = 'off';

                P1_sx_field.Visible = 'off';
                P1_s2_field.Visible = 'off';
                P2_sx_field.Visible = 'off';
                P2_s2_field.Visible = 'off';
                Px_sx_field.Visible = 'off';
                Px_s2_field.Visible = 'off';

                Px_T1_field.Enable = 'off';
                Px_nu1_field.Enable = 'off';
                Px_u1_field.Enable = 'off';
                Px_h1_field.Enable = 'off';
                Px_s1_field.Enable = 'off';

                Px_pressure.Enable = 'off';
            end
        elseif stepNumber == 3
            prop_type_panel.Visible = 'on';
            prop_type_panel.Enable = 'of';
            first_type_panel.Visible = 'on';
            first_type_panel.Enable = 'off';
            first_lookup_panel.Visible = 'on';
            first_lookup_panel.Enable = 'on';
            second_lookup_panel.Visible = 'off';

            if prop_nu.Value == 0
                first_lookup_rb_nu.Visible = 'off';
                second_lookup_rb_nu.Visible = 'off';
                if first_lookup_rb_nu.Value
                    first_lookup_rb_T.Value = 1;
                end
                if second_lookup_rb_nu.Value
                    second_lookup_rb_T.Value = 1;
                end
            else
                first_lookup_rb_nu.Visible = 'on';
            end

            if prop_u.Value == 0
                first_lookup_rb_u.Visible = 'off';
                second_lookup_rb_u.Visible = 'off';
                if first_lookup_rb_u.Value
                    first_lookup_rb_T.Value = 1;
                end
                if second_lookup_rb_u.Value
                    second_lookup_rb_T.Value = 1;
                end
            else
                first_lookup_rb_u.Visible = 'on';
            end

            if prop_h.Value == 0
                first_lookup_rb_h.Visible = 'off';
                second_lookup_rb_h.Visible = 'off';
                if first_lookup_rb_h.Value
                    first_lookup_rb_T.Value = 1;
                end
                if second_lookup_rb_h.Value
                    second_lookup_rb_T.Value = 1;
                end
            else
                first_lookup_rb_h.Visible = 'on';
            end

            if prop_s.Value == 0
                first_lookup_rb_s.Visible = 'off';
                second_lookup_rb_s.Visible = 'off';
                if first_lookup_rb_s.Value
                    first_lookup_rb_T.Value = 1;
                end
                if second_lookup_rb_s.Value
                    second_lookup_rb_T.Value = 1;
                end
            else
                first_lookup_rb_s.Visible = 'on';
            end

            if first_lookup_rb_T.Value
                P1_Tx_field.Enable = 'on';
                Px_T1_field.Enable = 'on';
                P1_nux_field.Enable = 'off';
                Px_nu1_field.Enable = 'off';
                P1_ux_field.Enable = 'off';
                Px_u1_field.Enable = 'off';
                P1_hx_field.Enable = 'off';
                Px_h1_field.Enable = 'off';
                P1_sx_field.Enable = 'off';
                Px_s1_field.Enable = 'off';
                Px_pressure.Enable = 'off';
            elseif first_lookup_rb_nu.Value
                P1_Tx_field.Enable = 'off';
                Px_T1_field.Enable = 'off';
                P1_nux_field.Enable = 'on';
                Px_nu1_field.Enable = 'on';
                P1_ux_field.Enable = 'off';
                Px_u1_field.Enable = 'off';
                P1_hx_field.Enable = 'off';
                Px_h1_field.Enable = 'off';
                P1_sx_field.Enable = 'off';
                Px_s1_field.Enable = 'off';
                Px_pressure.Enable = 'off';
            elseif first_lookup_rb_u.Value
                P1_Tx_field.Enable = 'off';
                Px_T1_field.Enable = 'off';
                P1_nux_field.Enable = 'off';
                Px_nu1_field.Enable = 'off';
                P1_ux_field.Enable = 'on';
                Px_u1_field.Enable = 'on';
                P1_hx_field.Enable = 'off';
                Px_h1_field.Enable = 'off';
                P1_sx_field.Enable = 'off';
                Px_s1_field.Enable = 'off';
                Px_pressure.Enable = 'off';
            elseif first_lookup_rb_h.Value
                P1_Tx_field.Enable = 'off';
                Px_T1_field.Enable = 'off';
                P1_nux_field.Enable = 'off';
                Px_nu1_field.Enable = 'off';
                P1_ux_field.Enable = 'off';
                Px_u1_field.Enable = 'off';
                P1_hx_field.Enable = 'on';
                Px_h1_field.Enable = 'on';
                P1_sx_field.Enable = 'off';
                Px_s1_field.Enable = 'off';
                Px_pressure.Enable = 'off';
            elseif first_lookup_rb_s.Value
                P1_Tx_field.Enable = 'off';
                Px_T1_field.Enable = 'off';
                P1_nux_field.Enable = 'off';
                Px_nu1_field.Enable = 'off';
                P1_ux_field.Enable = 'off';
                Px_u1_field.Enable = 'off';
                P1_hx_field.Enable = 'off';
                Px_h1_field.Enable = 'off';
                P1_sx_field.Enable = 'on';
                Px_s1_field.Enable = 'on';
                Px_pressure.Enable = 'off';
            else
                P1_Tx_field.Enable = 'off';
                Px_T1_field.Enable = 'off';
                P1_nux_field.Enable = 'off';
                Px_nu1_field.Enable = 'off';
                P1_ux_field.Enable = 'off';
                Px_u1_field.Enable = 'off';
                P1_hx_field.Enable = 'off';
                Px_h1_field.Enable = 'off';
                P1_sx_field.Enable = 'off';
                Px_s1_field.Enable = 'off';
                Px_pressure.Enable = 'on';
            end

            if first_type_rb_Pressure.Value
                if first_lookup_rb_nu.Value
                    Px_pressure.Value = interpolate(P1_nu1_field.Value, ...
                        P2_nu1_field.Value, P1_pressure.Value, ...
                        P2_pressure.Value, Px_nu1_field.Value);

                    Px_T1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_T1_field.Value, ...
                        P2_T1_field.Value, Px_pressure.Value);

                    Px_u1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_u1_field.Value, ...
                        P2_u1_field.Value, Px_pressure.Value);

                    Px_h1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_h1_field.Value, ...
                        P2_h1_field.Value, Px_pressure.Value);

                    Px_s1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_s1_field.Value, ...
                        P2_s1_field.Value, Px_pressure.Value);

                elseif first_lookup_rb_u.Value
                    Px_pressure.Value = interpolate(P1_u1_field.Value, ...
                        P2_u1_field.Value, P1_pressure.Value, ...
                        P2_pressure.Value, Px_u1_field.Value);

                    Px_T1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_T1_field.Value, ...
                        P2_T1_field.Value, Px_pressure.Value);

                    Px_nu1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_nu1_field.Value, ...
                        P2_nu1_field.Value, Px_pressure.Value);

                    Px_h1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_h1_field.Value, ...
                        P2_h1_field.Value, Px_pressure.Value);

                    Px_s1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_s1_field.Value, ...
                        P2_s1_field.Value, Px_pressure.Value);

                elseif first_lookup_rb_h.Value
                    Px_pressure.Value = interpolate(P1_h1_field.Value, ...
                        P2_h1_field.Value, P1_pressure.Value, ...
                        P2_pressure.Value, Px_h1_field.Value);

                    Px_T1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_T1_field.Value, ...
                        P2_T1_field.Value, Px_pressure.Value);

                    Px_nu1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_nu1_field.Value, ...
                        P2_nu1_field.Value, Px_pressure.Value);

                    Px_u1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_u1_field.Value, ...
                        P2_u1_field.Value, Px_pressure.Value);

                    Px_s1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_s1_field.Value, ...
                        P2_s1_field.Value, Px_pressure.Value);

                elseif first_lookup_rb_s.Value
                    Px_pressure.Value = interpolate(P1_s1_field.Value, ...
                        P2_s1_field.Value, P1_pressure.Value, ...
                        P2_pressure.Value, Px_s1_field.Value);

                    Px_T1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_T1_field.Value, ...
                        P2_T1_field.Value, Px_pressure.Value);

                    Px_nu1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_nu1_field.Value, ...
                        P2_nu1_field.Value, Px_pressure.Value);

                    Px_u1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_u1_field.Value, ...
                        P2_u1_field.Value, Px_pressure.Value);

                    Px_h1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_h1_field.Value, ...
                        P2_h1_field.Value, Px_pressure.Value);
                else
                    Px_T1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_T1_field.Value, ...
                        P2_T1_field.Value, Px_pressure.Value);

                    Px_nu1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_nu1_field.Value, ...
                        P2_nu1_field.Value, Px_pressure.Value);

                    Px_u1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_u1_field.Value, ...
                        P2_u1_field.Value, Px_pressure.Value);

                    Px_h1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_h1_field.Value, ...
                        P2_h1_field.Value, Px_pressure.Value);

                    Px_s1_field.Value = interpolate(P1_pressure.Value, ...
                        P2_pressure.Value, P1_s1_field.Value, ...
                        P2_s1_field.Value, Px_pressure.Value);
                end
            else
                if first_lookup_rb_T.Value
                    P1_nux_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_nu1_field.Value, ...
                        P1_nu2_field.Value, P1_Tx_field.Value);
                    P1_ux_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_u1_field.Value, ...
                        P1_u2_field.Value, P1_Tx_field.Value);
                    P1_hx_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_h1_field.Value, ...
                        P1_h2_field.Value, P1_Tx_field.Value);
                    P1_sx_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_s1_field.Value, ...
                        P1_s2_field.Value, P1_Tx_field.Value);
                elseif first_lookup_rb_nu.Value
                    P1_Tx_field.Value = interpolate(P1_nu1_field.Value, ...
                        P1_nu2_field.Value, P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_nux_field.Value);

                    P1_ux_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_u1_field.Value, ...
                        P1_u2_field.Value, P1_Tx_field.Value);
                    P1_hx_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_h1_field.Value, ...
                        P1_h2_field.Value, P1_Tx_field.Value);
                    P1_sx_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_s1_field.Value, ...
                        P1_s2_field.Value, P1_Tx_field.Value);
                elseif first_lookup_rb_u.Value
                    P1_Tx_field.Value = interpolate(P1_u1_field.Value, ...
                        P1_u2_field.Value, P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_ux_field.Value);

                    P1_nux_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_nu1_field.Value, ...
                        P1_nu2_field.Value, P1_Tx_field.Value);
                    P1_hx_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_h1_field.Value, ...
                        P1_h2_field.Value, P1_Tx_field.Value);
                    P1_sx_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_s1_field.Value, ...
                        P1_s2_field.Value, P1_Tx_field.Value);
                elseif first_lookup_rb_h.Value
                    P1_Tx_field.Value = interpolate(P1_h1_field.Value, ...
                        P1_h2_field.Value, P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_hx_field.Value);

                    P1_nux_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_nu1_field.Value, ...
                        P1_nu2_field.Value, P1_Tx_field.Value);
                    P1_ux_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_u1_field.Value, ...
                        P1_u2_field.Value, P1_Tx_field.Value);
                    P1_sx_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_s1_field.Value, ...
                        P1_s2_field.Value, P1_Tx_field.Value);
                else
                    P1_Tx_field.Value = interpolate(P1_s1_field.Value, ...
                        P1_s2_field.Value, P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_sx_field.Value);

                    P1_nux_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_nu1_field.Value, ...
                        P1_nu2_field.Value, P1_Tx_field.Value);
                    P1_ux_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_u1_field.Value, ...
                        P1_u2_field.Value, P1_Tx_field.Value);
                    P1_hx_field.Value = interpolate(P1_T1_field.Value, ...
                        P1_T2_field.Value, P1_h1_field.Value, ...
                        P1_h2_field.Value, P1_Tx_field.Value);
                end
            end

            if first_type_rb_other.Value
                first_lookup_rb_P.Visible = 'off';
                first_lookup_rb_T.Visible = 'on';
                if first_lookup_rb_P.Value
                    first_lookup_rb_T.Value = 1;
                end
                P2_panel.Visible = 'off';
                Px_panel.Visible = 'off';




                if (prop_nu.Value + prop_u.Value + prop_h.Value + prop_s.Value) == 0
                    prop_h.Value = 1;
                end

                first_lookup_rb_h.Visible = 'on';
                first_lookup_rb_nu.Visible = 'on';
                first_lookup_rb_u.Visible = 'on';
                first_lookup_rb_h.Visible = 'on';
                first_lookup_rb_s.Visible = 'on';

                second_lookup_rb_h.Visible = 'on';
                second_lookup_rb_nu.Visible = 'on';
                second_lookup_rb_u.Visible = 'on';
                second_lookup_rb_h.Visible = 'on';
                second_lookup_rb_s.Visible = 'on';

                if prop_nu.Value == 0
                    first_lookup_rb_nu.Visible = 'off';
                    second_lookup_rb_nu.Visible = 'off';

                    P1_nu1_field.Visible = 'off';
                    P1_nux_field.Visible = 'off';
                    P1_nu2_field.Visible = 'off';
                    P2_nu1_field.Visible = 'off';
                    P2_nux_field.Visible = 'off';
                    P2_nu2_field.Visible = 'off';
                    Px_nu1_field.Visible = 'off';
                    Px_nux_field.Visible = 'off';
                    Px_nu2_field.Visible = 'off';
                else
                    P1_nu1_field.Visible = 'on';
                    P1_nux_field.Visible = 'on';
                    P1_nu2_field.Visible = 'on';
                    P2_nu1_field.Visible = 'on';
                    P2_nux_field.Visible = 'on';
                    P2_nu2_field.Visible = 'on';
                    Px_nu1_field.Visible = 'on';
                    Px_nux_field.Visible = 'on';
                    Px_nu2_field.Visible = 'on';
                end

                if prop_u.Value == 0
                    first_lookup_rb_u.Visible = 'off';
                    second_lookup_rb_u.Visible = 'off';

                    P1_u1_field.Visible = 'off';
                    P1_ux_field.Visible = 'off';
                    P1_u2_field.Visible = 'off';
                    P2_u1_field.Visible = 'off';
                    P2_ux_field.Visible = 'off';
                    P2_u2_field.Visible = 'off';
                    Px_u1_field.Visible = 'off';
                    Px_ux_field.Visible = 'off';
                    Px_u2_field.Visible = 'off';
                else
                    P1_u1_field.Visible = 'on';
                    P1_ux_field.Visible = 'on';
                    P1_u2_field.Visible = 'on';
                    P2_u1_field.Visible = 'on';
                    P2_ux_field.Visible = 'on';
                    P2_u2_field.Visible = 'on';
                    Px_u1_field.Visible = 'on';
                    Px_ux_field.Visible = 'on';
                    Px_u2_field.Visible = 'on';
                end

                if prop_h.Value == 0

                    first_lookup_rb_h.Visible = 'off';
                    second_lookup_rb_h.Visible = 'off';

                    P1_h1_field.Visible = 'off';
                    P1_hx_field.Visible = 'off';
                    P1_h2_field.Visible = 'off';
                    P2_h1_field.Visible = 'off';
                    P2_hx_field.Visible = 'off';
                    P2_h2_field.Visible = 'off';
                    Px_h1_field.Visible = 'off';
                    Px_hx_field.Visible = 'off';
                    Px_h2_field.Visible = 'off';
                else
                    P1_h1_field.Visible = 'on';
                    P1_hx_field.Visible = 'on';
                    P1_h2_field.Visible = 'on';
                    P2_h1_field.Visible = 'on';
                    P2_hx_field.Visible = 'on';
                    P2_h2_field.Visible = 'on';
                    Px_h1_field.Visible = 'on';
                    Px_hx_field.Visible = 'on';
                    Px_h2_field.Visible = 'on';
                end

                if prop_s.Value == 0

                    first_lookup_rb_s.Visible = 'off';
                    second_lookup_rb_s.Visible = 'off';

                    P1_s1_field.Visible = 'off';
                    P1_sx_field.Visible = 'off';
                    P1_s2_field.Visible = 'off';
                    P2_s1_field.Visible = 'off';
                    P2_sx_field.Visible = 'off';
                    P2_s2_field.Visible = 'off';
                    Px_s1_field.Visible = 'off';
                    Px_sx_field.Visible = 'off';
                    Px_s2_field.Visible = 'off';
                else
                    P1_s1_field.Visible = 'on';
                    P1_sx_field.Visible = 'on';
                    P1_s2_field.Visible = 'on';
                    P2_s1_field.Visible = 'on';
                    P2_sx_field.Visible = 'on';
                    P2_s2_field.Visible = 'on';
                    Px_s1_field.Visible = 'on';
                    Px_sx_field.Visible = 'on';
                    Px_s2_field.Visible = 'on';
                end
                P1_Tx_field.Visible = 'on';
                P1_T2_field.Visible = 'on';
                P2_Tx_field.Visible = 'on';
                P2_T2_field.Visible = 'on';
                Px_Tx_field.Visible = 'on';
                Px_T2_field.Visible = 'on';

            else
                first_lookup_rb_T.Visible = 'off';
                first_lookup_rb_P.Visible = 'on';
                if first_lookup_rb_T.Value
                    first_lookup_rb_P.Value = 1;
                end

                P2_panel.Visible = 'on';
                Px_panel.Visible = 'on';

                P1_Tx_field.Visible = 'off';
                P1_T2_field.Visible = 'off';
                P2_Tx_field.Visible = 'off';
                P2_T2_field.Visible = 'off';
                Px_Tx_field.Visible = 'off';
                Px_T2_field.Visible = 'off';

                P1_nux_field.Visible = 'off';
                P1_nu2_field.Visible = 'off';
                P2_nux_field.Visible = 'off';
                P2_nu2_field.Visible = 'off';
                Px_nux_field.Visible = 'off';
                Px_nu2_field.Visible = 'off';

                P1_ux_field.Visible = 'off';
                P1_u2_field.Visible = 'off';
                P2_ux_field.Visible = 'off';
                P2_u2_field.Visible = 'off';
                Px_ux_field.Visible = 'off';
                Px_u2_field.Visible = 'off';

                P1_hx_field.Visible = 'off';
                P1_h2_field.Visible = 'off';
                P2_hx_field.Visible = 'off';
                P2_h2_field.Visible = 'off';
                Px_hx_field.Visible = 'off';
                Px_h2_field.Visible = 'off';

                P1_sx_field.Visible = 'off';
                P1_s2_field.Visible = 'off';
                P2_sx_field.Visible = 'off';
                P2_s2_field.Visible = 'off';
                Px_sx_field.Visible = 'off';
                Px_s2_field.Visible = 'off';

            end

        else % STEP 4

            % Change visability of input panels
            prop_type_panel.Visible = 'on';
            prop_type_panel.Enable = 'of';
            first_type_panel.Visible = 'on';
            first_type_panel.Enable = 'off';
            first_lookup_panel.Visible = 'on';
            first_lookup_panel.Enable = 'off';
            second_lookup_panel.Visible = 'on';
            second_lookup_panel.Enable = 'on';

            % Make all pressure panels visible
            Px_panel.Visible = 'on';
            P2_panel.Visible = 'on';

            % Turn all previous enabled interpolation values off
            Px_T1_field.Enable = 'off';
            Px_nu1_field.Enable = 'off';
            Px_u1_field.Enable = 'off';
            Px_h1_field.Enable = 'off';
            Px_s1_field.Enable = 'off';

            P1_Tx_field.Enable = 'off';
            P1_nux_field.Enable = 'off';
            P1_ux_field.Enable = 'off';
            P1_hx_field.Enable = 'off';
            P1_sx_field.Enable = 'off';

            % Turn off rest of Px panel
            Px_Tx_field.Enable = 'off';
            Px_nux_field.Enable = 'off';
            Px_ux_field.Enable = 'off';
            Px_hx_field.Enable = 'off';
            Px_sx_field.Enable = 'off';

            Px_T2_field.Enable = 'off';
            Px_nu2_field.Enable = 'off';
            Px_u2_field.Enable = 'off';
            Px_h2_field.Enable = 'off';
            Px_s2_field.Enable = 'off';

            Px_pressure.Enable = 'off';

            % Turn all field visibilities on
            P1_Tx_field.Visible = 'on';
            P1_nux_field.Visible = 'on';
            P1_ux_field.Visible = 'on';
            P1_hx_field.Visible = 'on';
            P1_sx_field.Visible = 'on';

            P1_T2_field.Visible = 'on';
            P1_nu2_field.Visible = 'on';
            P1_u2_field.Visible = 'on';
            P1_h2_field.Visible = 'on';
            P1_s2_field.Visible = 'on';

            P2_Tx_field.Visible = 'on';
            P2_nux_field.Visible = 'on';
            P2_ux_field.Visible = 'on';
            P2_hx_field.Visible = 'on';
            P2_sx_field.Visible = 'on';

            P2_T2_field.Visible = 'on';
            P2_nu2_field.Visible = 'on';
            P2_u2_field.Visible = 'on';
            P2_h2_field.Visible = 'on';
            P2_s2_field.Visible = 'on';

            Px_Tx_field.Visible = 'on';
            Px_nux_field.Visible = 'on';
            Px_ux_field.Visible = 'on';
            Px_hx_field.Visible = 'on';
            Px_sx_field.Visible = 'on';

            Px_T2_field.Visible = 'on';
            Px_nu2_field.Visible = 'on';
            Px_u2_field.Visible = 'on';
            Px_h2_field.Visible = 'on';
            Px_s2_field.Visible = 'on';

            % Control visibility of radio buttons

            % Turn on all buttons that should be
            second_lookup_rb_T.Visible = 'on';
            second_lookup_rb_P.Visible = 'on';
            if prop_nu.Value
                second_lookup_rb_nu.Visible = 'on';
            end
            if prop_u.Value
                second_lookup_rb_u.Visible = 'on';
            end
            if prop_h.Value
                second_lookup_rb_h.Visible = 'on';
            end
            if prop_s.Value
                second_lookup_rb_s.Visible = 'on';
            end

            % Turn off any button that was selected last time
            if first_lookup_rb_T.Value % If first lookup was temperature
                second_lookup_rb_T.Visible = 'off'; % Turn off temperature option the second time
                % Switch radio button to visible, available option if
                % change is needed
                if second_lookup_rb_T.Value
                    if prop_nu.Value
                        second_lookup_rb_nu.Value = 1;
                    elseif prop_u.Value
                        second_lookup_rb_u.Value = 1;
                    elseif prop_h.Value
                        second_lookup_rb_h.Value = 1;
                    elseif prop_s.Value
                        second_lookup_rb_s.Value = 1;
                    else
                        second_lookup_rb_P.Value = 1;
                    end
                end
            elseif first_lookup_rb_nu.Value % If first lookup was specific volume
                second_lookup_rb_nu.Visible = 'off'; % Turn off specific volume option the second time
                % Switch radio button to temperature if change is needed
                if second_lookup_rb_nu.Value
                    second_lookup_rb_T.Value = 1;
                end
            elseif first_lookup_rb_u.Value % If first lookup was internal energy
                second_lookup_rb_u.Visible = 'off'; % Turn off internal option the second time
                % Switch radio button to temperature if change is needed
                if second_lookup_rb_u.Value
                    second_lookup_rb_T.Value = 1;
                end
            elseif first_lookup_rb_h.Value % If first lookup was enthalpy
                second_lookup_rb_h.Visible = 'off'; % Turn off enthalpy option the second time
                % Switch radio button to temperature if change is needed
                if second_lookup_rb_h.Value
                    second_lookup_rb_T.Value = 1;
                end
            elseif first_lookup_rb_s.Value % If first lookup was entropy
                second_lookup_rb_s.Visible = 'off'; % Turn off entropy option the second time
                % Switch radio button to temperature if change is needed
                if second_lookup_rb_s.Value
                    second_lookup_rb_T.Value = 1;
                end
            else % If first lookup was pressure
                second_lookup_rb_P.Visible = 'off'; % Turn off entropy option the second time
                % Switch radio button to temperature if change is needed
                if second_lookup_rb_P.Value
                    second_lookup_rb_T.Value = 1;
                end
            end

            if first_type_rb_other.Value % If it was constant pressure
                if first_lookup_rb_T.Value % If the first lookup was temperature
                    P2_Tx_field.Value = P1_Tx_field.Value; % Set second panel Tx to input value

                    % Interpolate the second panel based on temperature
                    P2_nux_field.Value = interpolate(P2_T1_field.Value, ...
                        P2_T2_field.Value, P2_nu1_field.Value, ...
                        P2_nu2_field.Value, P2_Tx_field.Value);
                    P2_ux_field.Value = interpolate(P2_T1_field.Value, ...
                        P2_T2_field.Value, P2_u1_field.Value, ...
                        P2_u2_field.Value, P2_Tx_field.Value);
                    P2_hx_field.Value = interpolate(P2_T1_field.Value, ...
                        P2_T2_field.Value, P2_h1_field.Value, ...
                        P2_h2_field.Value, P2_Tx_field.Value);
                    P2_sx_field.Value = interpolate(P2_T1_field.Value, ...
                        P2_T2_field.Value, P2_s1_field.Value, ...
                        P2_s2_field.Value, P2_Tx_field.Value);
                elseif first_lookup_rb_nu.Value % If the first lookup was specific volume
                    P2_nux_field.Value = P1_nux_field.Value; % Set second panel nux to input value

                    % Interpolate the second panel based on specific volume
                    P2_Tx_field.Value = interpolate(P2_nu1_field.Value, ...
                        P2_nu2_field.Value, P2_T1_field.Value, ...
                        P2_T2_field.Value, P2_nux_field.Value);
                    P2_ux_field.Value = interpolate(P2_nu1_field.Value, ...
                        P2_nu2_field.Value, P2_u1_field.Value, ...
                        P2_u2_field.Value, P2_nux_field.Value);
                    P2_hx_field.Value = interpolate(P2_nu1_field.Value, ...
                        P2_nu2_field.Value, P2_h1_field.Value, ...
                        P2_h2_field.Value, P2_nux_field.Value);
                    P2_sx_field.Value = interpolate(P2_nu1_field.Value, ...
                        P2_nu2_field.Value, P2_s1_field.Value, ...
                        P2_s2_field.Value, P2_nux_field.Value);
                elseif first_lookup_rb_u.Value % If the first lookup was internal energy
                    P2_ux_field.Value = P1_ux_field.Value; % Set second panel ux to input value

                    % Interpolate the second panel based on internal energy
                    P2_Tx_field.Value = interpolate(P2_u1_field.Value, ...
                        P2_u2_field.Value, P2_T1_field.Value, ...
                        P2_T2_field.Value, P2_ux_field.Value);
                    P2_nux_field.Value = interpolate(P2_u1_field.Value, ...
                        P2_u2_field.Value, P2_nu1_field.Value, ...
                        P2_nu2_field.Value, P2_ux_field.Value);
                    P2_hx_field.Value = interpolate(P2_u1_field.Value, ...
                        P2_u2_field.Value, P2_h1_field.Value, ...
                        P2_h2_field.Value, P2_ux_field.Value);
                    P2_sx_field.Value = interpolate(P2_u1_field.Value, ...
                        P2_u2_field.Value, P2_s1_field.Value, ...
                        P2_s2_field.Value, P2_ux_field.Value);
                elseif first_lookup_rb_h.Value % If the first lookup was enthalpy
                    P2_hx_field.Value = P1_hx_field.Value; % Set second panel hx to input value

                    % Interpolate the second panel based on enthalpy
                    P2_Tx_field.Value = interpolate(P2_h1_field.Value, ...
                        P2_h2_field.Value, P2_T1_field.Value, ...
                        P2_T2_field.Value, P2_hx_field.Value);
                    P2_nux_field.Value = interpolate(P2_h1_field.Value, ...
                        P2_h2_field.Value, P2_nu1_field.Value, ...
                        P2_nu2_field.Value, P2_hx_field.Value);
                    P2_ux_field.Value = interpolate(P2_h1_field.Value, ...
                        P2_h2_field.Value, P2_u1_field.Value, ...
                        P2_u2_field.Value, P2_hx_field.Value);
                    P2_sx_field.Value = interpolate(P2_h1_field.Value, ...
                        P2_h2_field.Value, P2_s1_field.Value, ...
                        P2_s2_field.Value, P2_hx_field.Value);
                else % If the first lookup was entropy
                    P2_sx_field.Value = P1_sx_field.Value; % Set second panel sx to input value

                    % Interpolate the second panel based on entropy
                    P2_Tx_field.Value = interpolate(P2_s1_field.Value, ...
                        P2_s2_field.Value, P2_T1_field.Value, ...
                        P2_T2_field.Value, P2_sx_field.Value);
                    P2_nux_field.Value = interpolate(P2_s1_field.Value, ...
                        P2_s2_field.Value, P2_nu1_field.Value, ...
                        P2_nu2_field.Value, P2_sx_field.Value);
                    P2_ux_field.Value = interpolate(P2_s1_field.Value, ...
                        P2_s2_field.Value, P2_u1_field.Value, ...
                        P2_u2_field.Value, P2_sx_field.Value);
                    P2_hx_field.Value = interpolate(P2_s1_field.Value, ...
                        P2_s2_field.Value, P2_h1_field.Value, ...
                        P2_h2_field.Value, P2_sx_field.Value);
                end

                % Interpolate x panel based on second input
                if second_lookup_rb_P.Value % If the second lookup was pressure
                    Px_pressure.Enable = 'on';
                elseif second_lookup_rb_T.Value % If second lookup was temperature
                    Px_Tx_field.Enable = 'on';
                    % Interpolate pressure from temperature
                    Px_pressure.Value = interpolate(P1_Tx_field.Value, ...
                        P2_Tx_field.Value, P1_pressure.Value, ...
                        P2_pressure.Value, Px_Tx_field.Value);
                elseif second_lookup_rb_nu.Value % If the second lookup was specific volume
                    Px_nux_field.Enable = 'on';
                    % Interpolate pressure from nu
                    Px_pressure.Value = interpolate(P1_nux_field.Value, ...
                        P2_nux_field.Value, P1_pressure.Value, ...
                        P2_pressure.Value, Px_nux_field.Value);
                elseif second_lookup_rb_u.Value % If second lookup was internal energy
                    Px_ux_field.Enable = 'on';
                    % Interpolate pressure from u
                    Px_pressure.Value = interpolate(P1_ux_field.Value, ...
                        P2_ux_field.Value, P1_pressure.Value, ...
                        P2_pressure.Value, Px_ux_field.Value);
                elseif second_lookup_rb_h.Value % If second lookup was enthalpy
                    Px_hx_field.Enable = 'on';
                    % Interpolate pressure from h
                    Px_pressure.Value = interpolate(P1_hx_field.Value, ...
                        P2_hx_field.Value, P1_pressure.Value, ...
                        P2_pressure.Value, Px_hx_field.Value);
                else % If second lookup was entropy
                    Px_sx_field.Enable = 'on';
                    % Interpolate pressure from s
                    Px_pressure.Value = interpolate(P1_sx_field.Value, ...
                        P2_sx_field.Value, P1_pressure.Value, ...
                        P2_pressure.Value, Px_sx_field.Value);
                end

                % Interpolate x panel based on pressure
                Px_T1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_T1_field.Value, ...
                    P2_T1_field.Value, Px_pressure.Value);
                Px_nu1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_nu1_field.Value, ...
                    P2_nu1_field.Value, Px_pressure.Value);
                Px_u1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_u1_field.Value, ...
                    P2_u1_field.Value, Px_pressure.Value);
                Px_h1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_h1_field.Value, ...
                    P2_h1_field.Value, Px_pressure.Value);
                Px_s1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_s1_field.Value, ...
                    P2_s1_field.Value, Px_pressure.Value);

                Px_Tx_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_Tx_field.Value, ...
                    P2_Tx_field.Value, Px_pressure.Value);
                Px_nux_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_nux_field.Value, ...
                    P2_nux_field.Value, Px_pressure.Value);
                Px_ux_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_ux_field.Value, ...
                    P2_ux_field.Value, Px_pressure.Value);
                Px_hx_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_hx_field.Value, ...
                    P2_hx_field.Value, Px_pressure.Value);
                Px_sx_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_sx_field.Value, ...
                    P2_sx_field.Value, Px_pressure.Value);

                Px_T2_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_T2_field.Value, ...
                    P2_T2_field.Value, Px_pressure.Value);
                Px_nu2_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_nu2_field.Value, ...
                    P2_nu2_field.Value, Px_pressure.Value);
                Px_u2_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_u2_field.Value, ...
                    P2_u2_field.Value, Px_pressure.Value);
                Px_h2_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_h2_field.Value, ...
                    P2_h2_field.Value, Px_pressure.Value);
                Px_s2_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_s2_field.Value, ...
                    P2_s2_field.Value, Px_pressure.Value);
            else % If it was cross-pressure
                % Interpolate Px values 1 and 2 (in case properties were
                % changed in this step)
                Px_T1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_T1_field.Value, ...
                    P2_T1_field.Value, Px_pressure.Value);
                Px_nu1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_nu1_field.Value, ...
                    P2_nu1_field.Value, Px_pressure.Value);
                Px_u1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_u1_field.Value, ...
                    P2_u1_field.Value, Px_pressure.Value);
                Px_h1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_h1_field.Value, ...
                    P2_h1_field.Value, Px_pressure.Value);
                Px_s1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_s1_field.Value, ...
                    P2_s1_field.Value, Px_pressure.Value);

                Px_T2_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_T2_field.Value, ...
                    P2_T2_field.Value, Px_pressure.Value);
                Px_nu1_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_nu2_field.Value, ...
                    P2_nu2_field.Value, Px_pressure.Value);
                Px_u2_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_u2_field.Value, ...
                    P2_u2_field.Value, Px_pressure.Value);
                Px_h2_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_h2_field.Value, ...
                    P2_h2_field.Value, Px_pressure.Value);
                Px_s2_field.Value = interpolate(P1_pressure.Value, ...
                    P2_pressure.Value, P1_s2_field.Value, ...
                    P2_s2_field.Value, Px_pressure.Value);
                if second_lookup_rb_T.Value % If second lookup is temperature
                    Px_Tx_field.Enable = 'on';
                elseif second_lookup_rb_nu.Value % If second lookup was specific volume
                    Px_nux_field.Enable = 'on';

                    % Interpolate temperature value for panel x based on
                    % input property
                    Px_Tx_field.Value = interpolate(Px_nu1_field.Value, ...
                        Px_nu2_field.Value, Px_T1_field.Value, ...
                        Px_T2_field.Value, Px_nux_field.Value);
                elseif second_lookup_rb_u.Value % If second lookup was internal energy
                    Px_ux_field.Enable = 'on';

                    % Interpolate temperature value for panel x based on
                    % input property
                    Px_Tx_field.Value = interpolate(Px_u1_field.Value, ...
                        Px_u2_field.Value, Px_T1_field.Value, ...
                        Px_T2_field.Value, Px_ux_field.Value);
                elseif second_lookup_rb_h.Value % If second lookup was enthalpy
                    Px_hx_field.Enable = 'on';

                    % Interpolate temperature value for panel x based on
                    % input property
                    Px_Tx_field.Value = interpolate(Px_h1_field.Value, ...
                        Px_h2_field.Value, Px_T1_field.Value, ...
                        Px_T2_field.Value, Px_hx_field.Value);
                else % If second lookup was entropy
                    Px_sx_field.Enable = 'on';

                    % Interpolate temperature value for panel x based on
                    % input property
                    Px_Tx_field.Value = interpolate(Px_s1_field.Value, ...
                        Px_s2_field.Value, Px_T1_field.Value, ...
                        Px_T2_field.Value, Px_sx_field.Value);
                end

                % Set first and second panels to have same temperature
                % value as panel x
                P1_Tx_field.Value = Px_Tx_field.Value;
                P2_Tx_field.Value = Px_Tx_field.Value;

                % Interpolate all panels based on temperature
                P1_nux_field.Value = interpolate(P1_T1_field.Value, ...
                    P1_T2_field.Value, P1_nu1_field.Value, ...
                    P1_nu2_field.Value, P1_Tx_field.Value);
                P1_ux_field.Value = interpolate(P1_T1_field.Value, ...
                    P1_T2_field.Value, P1_u1_field.Value, ...
                    P1_u2_field.Value, P1_Tx_field.Value);
                P1_hx_field.Value = interpolate(P1_T1_field.Value, ...
                    P1_T2_field.Value, P1_h1_field.Value, ...
                    P1_h2_field.Value, P1_Tx_field.Value);
                P1_sx_field.Value = interpolate(P1_T1_field.Value, ...
                    P1_T2_field.Value, P1_s1_field.Value, ...
                    P1_s2_field.Value, P1_Tx_field.Value);

                P2_nux_field.Value = interpolate(P2_T1_field.Value, ...
                    P2_T2_field.Value, P2_nu1_field.Value, ...
                    P2_nu2_field.Value, P2_Tx_field.Value);
                P2_ux_field.Value = interpolate(P2_T1_field.Value, ...
                    P2_T2_field.Value, P2_u1_field.Value, ...
                    P2_u2_field.Value, P2_Tx_field.Value);
                P2_hx_field.Value = interpolate(P2_T1_field.Value, ...
                    P2_T2_field.Value, P2_h1_field.Value, ...
                    P2_h2_field.Value, P2_Tx_field.Value);
                P2_sx_field.Value = interpolate(P2_T1_field.Value, ...
                    P2_T2_field.Value, P2_s1_field.Value, ...
                    P2_s2_field.Value, P2_Tx_field.Value);

                Px_nux_field.Value = interpolate(Px_T1_field.Value, ...
                    Px_T2_field.Value, Px_nu1_field.Value, ...
                    Px_nu2_field.Value, Px_Tx_field.Value);
                Px_ux_field.Value = interpolate(Px_T1_field.Value, ...
                    Px_T2_field.Value, Px_u1_field.Value, ...
                    Px_u2_field.Value, Px_Tx_field.Value);
                Px_hx_field.Value = interpolate(Px_T1_field.Value, ...
                    Px_T2_field.Value, Px_h1_field.Value, ...
                    Px_h2_field.Value, Px_Tx_field.Value);
                Px_sx_field.Value = interpolate(Px_T1_field.Value, ...
                    Px_T2_field.Value, Px_s1_field.Value, ...
                    Px_s2_field.Value, Px_Tx_field.Value);
            end

            % Turn off nonexistant property fields
            if prop_nu.Value == 0
                first_lookup_rb_nu.Visible = 'off';
                second_lookup_rb_nu.Visible = 'off';

                P1_nu1_field.Visible = 'off';
                P1_nux_field.Visible = 'off';
                P1_nu2_field.Visible = 'off';
                P2_nu1_field.Visible = 'off';
                P2_nux_field.Visible = 'off';
                P2_nu2_field.Visible = 'off';
                Px_nu1_field.Visible = 'off';
                Px_nux_field.Visible = 'off';
                Px_nu2_field.Visible = 'off';
            else
                P1_nu1_field.Visible = 'on';
                P1_nux_field.Visible = 'on';
                P1_nu2_field.Visible = 'on';
                P2_nu1_field.Visible = 'on';
                P2_nux_field.Visible = 'on';
                P2_nu2_field.Visible = 'on';
                Px_nu1_field.Visible = 'on';
                Px_nux_field.Visible = 'on';
                Px_nu2_field.Visible = 'on';
            end

            if prop_u.Value == 0

                first_lookup_rb_u.Visible = 'off';
                second_lookup_rb_u.Visible = 'off';

                P1_u1_field.Visible = 'off';
                P1_ux_field.Visible = 'off';
                P1_u2_field.Visible = 'off';
                P2_u1_field.Visible = 'off';
                P2_ux_field.Visible = 'off';
                P2_u2_field.Visible = 'off';
                Px_u1_field.Visible = 'off';
                Px_ux_field.Visible = 'off';
                Px_u2_field.Visible = 'off';
            else
                P1_u1_field.Visible = 'on';
                P1_ux_field.Visible = 'on';
                P1_u2_field.Visible = 'on';
                P2_u1_field.Visible = 'on';
                P2_ux_field.Visible = 'on';
                P2_u2_field.Visible = 'on';
                Px_u1_field.Visible = 'on';
                Px_ux_field.Visible = 'on';
                Px_u2_field.Visible = 'on';
            end

            if prop_h.Value == 0

                first_lookup_rb_h.Visible = 'off';
                second_lookup_rb_h.Visible = 'off';

                P1_h1_field.Visible = 'off';
                P1_hx_field.Visible = 'off';
                P1_h2_field.Visible = 'off';
                P2_h1_field.Visible = 'off';
                P2_hx_field.Visible = 'off';
                P2_h2_field.Visible = 'off';
                Px_h1_field.Visible = 'off';
                Px_hx_field.Visible = 'off';
                Px_h2_field.Visible = 'off';
            else
                P1_h1_field.Visible = 'on';
                P1_hx_field.Visible = 'on';
                P1_h2_field.Visible = 'on';
                P2_h1_field.Visible = 'on';
                P2_hx_field.Visible = 'on';
                P2_h2_field.Visible = 'on';
                Px_h1_field.Visible = 'on';
                Px_hx_field.Visible = 'on';
                Px_h2_field.Visible = 'on';
            end

            if prop_s.Value == 0

                first_lookup_rb_s.Visible = 'off';
                second_lookup_rb_s.Visible = 'off';

                P1_s1_field.Visible = 'off';
                P1_sx_field.Visible = 'off';
                P1_s2_field.Visible = 'off';
                P2_s1_field.Visible = 'off';
                P2_sx_field.Visible = 'off';
                P2_s2_field.Visible = 'off';
                Px_s1_field.Visible = 'off';
                Px_sx_field.Visible = 'off';
                Px_s2_field.Visible = 'off';
            else
                P1_s1_field.Visible = 'on';
                P1_sx_field.Visible = 'on';
                P1_s2_field.Visible = 'on';
                P2_s1_field.Visible = 'on';
                P2_sx_field.Visible = 'on';
                P2_s2_field.Visible = 'on';
                Px_s1_field.Visible = 'on';
                Px_sx_field.Visible = 'on';
                Px_s2_field.Visible = 'on';
            end
        end
    end
end