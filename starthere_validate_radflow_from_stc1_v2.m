%%
clc;

[resA, resB, data] = validate_radflow_from_stc1_v2('C:\zy\code\mcode\adofflineanalysis-main\MSTd.mat','sinkin','ves');

% Plot Nearest-4 angles:
plot_radflow_agreement(data.theta_deg, data.Rfull_z, resA.ang_w, 'sinkin', 'MSTd VIS — Nearest-4');

% Plot Fit-then-fold angles:
plot_radflow_agreement(data.theta_deg, data.Rfull_z, resB.ang_w, 'sinkin', 'MSTd VIS — Fit-then-fold');


%%

% Visual condition, sink-in geometry
[res_vis, data_vis] = validate_RSI_from_stc1('E:\zy\mcode\adofflineanalysis-main\MSTd.mat','sinkin','vis', true);

% Vestibular condition, center-out geometry
[res_ves, data_ves] = validate_RSI_from_stc1('E:\zy\mcode\adofflineanalysis-main\MSTd.mat','centerout','ves', true);
