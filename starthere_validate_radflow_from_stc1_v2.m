%%
clc;

[resA, resB, data] = validate_radflow_from_stc1_v2('C:\zy\code\mcode\adofflineanalysis-main\MSTd.mat','sinkin','ves');

% Plot Nearest-4 angles:
plot_radflow_agreement(data.theta_deg, data.Rfull_z, resA.ang_w, 'sinkin', 'MSTd VIS — Nearest-4');

% Plot Fit-then-fold angles:
plot_radflow_agreement(data.theta_deg, data.Rfull_z, resB.ang_w, 'sinkin', 'MSTd VIS — Fit-then-fold');


