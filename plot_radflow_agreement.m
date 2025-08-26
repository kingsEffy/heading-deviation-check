function plot_radflow_agreement(theta_deg, Rfull, ang_w, stim_mode, figtitle)
% Visualize accuracy: dataset "reported" headings vs headings from your script.
%
% Inputs
%   theta_deg : 1×K stimulus headings (deg)
%   Rfull     : nUnits×K trial-averaged responses for those headings
%   ang_w     : nUnits×1 angles from linearfit_radical (radians)
%   stim_mode : 'sinkin' or 'centerout'  (affects 180° alignment check)
%   figtitle  : (optional) title string
%
% What it draws
%   (A) wrap-aware scatter of reported vs computed angles (+ unity line)
%   (B) circular error rose (bins=10°)
%   (C) CDF of |error|  (read off % within thresholds)
%   (D) polar "spokes": reported → computed (subset to avoid clutter)
%
% Printed metrics
%   Circular MAE, median, % within ±22.5° & ±11.25°, circular corr (Fisher–Lee ρ)

    if nargin < 5, figtitle = 'Heading agreement'; end
    % 1) "Reported" heading per unit = argmax across ALL headings (zscore for shape)
    Rz = zscore(Rfull, 0, 2);
    [~, ixMax] = max(Rz, [], 2);
    head_rep = wrap360(theta_deg(ixMax));

    % 2) Your computed headings (deg), aligned to the dataset
    ang_deg = wrap360(rad2deg(ang_w));
    ang_aligned = align_angles_to_truth(ang_deg, head_rep, strcmpi(stim_mode,'centerout'));

    % 3) Circular errors
    err = wrap180(ang_aligned - head_rep);
    absErr = abs(err);

    % 4) Summary metrics
    mae   = mean(absErr);
    med   = median(absErr);
    pct22 = mean(absErr <= 22.5)*100;
    pct11 = mean(absErr <= 11.25)*100;
    rho   = circ_corr_deg(ang_aligned, head_rep);

    % ---------- Figure ----------
    f = figure('Color','w','Position',[80 80 1200 800]);
    tl = tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
    title(tl, sprintf('%s  |  MAE=%.2f°, med=%.2f° | ≤22.5°=%.1f%%, ≤11.25°=%.1f%% | ρ=%.3f', ...
                      figtitle, mae, med, pct22, pct11, rho), 'FontWeight','bold');

    % (A) Wrap-aware scatter (reported vs computed) + unity
    nexttile;
    [x,y] = wrapaware_xy(head_rep, ang_aligned);
    scatter(x, y, 16, 'filled', 'MarkerFaceAlpha', 0.7); hold on;
    plot([0 360],[0 360],'k--','LineWidth',1);
    xlim([0 360]); ylim([0 360]);
    xticks(0:45:360); yticks(0:45:360);
    xlabel('Reported heading (deg)'); ylabel('Computed heading (deg)');
    axis square; grid on; title('Reported vs computed (wrap-aware)');

    % (B) Circular error rose (10° bins)
    nexttile;
    edges = deg2rad(-180:10:180);
    polarhistogram(deg2rad(err), edges, 'Normalization','probability');
    title('Circular error (deg)'); pax = gca; pax.ThetaZeroLocation='top'; pax.ThetaDir='clockwise';
    % Annotate
    text(0.02, 0.98, sprintf('MAE=%.2f°  med=%.2f°', mae, med), 'Units','normalized', ...
         'HorizontalAlignment','left','VerticalAlignment','top');

    % (C) CDF of |error|
    nexttile;
    % [F,xx] = ecdf(absErr);
    % plot(xx, F, 'LineWidth', 2); hold on;
    % xline(22.5, ':'); xline(11.25, ':');
    % y22 = interp1(xx, F, 22.5, 'linear','extrap')*100;
    % y11 = interp1(xx, F, 11.25,'linear','extrap')*100;
    % grid on; xlim([0 180]); ylim([0 1]);
    % xlabel('|error| (deg)'); ylabel('CDF');
    % title(sprintf('CDF  |  ≤22.5°=%.1f%%, ≤11.25°=%.1f%%', y22, y11));
    % (C) CDF of |error|
    valid = ~isnan(absErr);
    [F,xx] = ecdf(absErr(valid));
    % ecdf is a step function; use stairs for clarity and no uniqueness issues
    stairs(xx, F, 'LineWidth', 2); hold on;
    
    % Exact percentages (no interp1 needed)
    y22 = mean(absErr(valid) <= 22.5)*100;
    y11 = mean(absErr(valid) <= 11.25)*100;
    
    % Visual guides and markers
    xline(22.5, ':'); xline(11.25, ':');
    plot(22.5, y22/100, 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
    plot(11.25, y11/100, 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
    grid on; xlim([0 180]); ylim([0 1]);
    xlabel('|error| (deg)'); ylabel('CDF');
    title(sprintf('CDF  |  ≤22.5°=%.1f%%, ≤11.25°=%.1f%%', y22, y11));

    % (D) Polar "spokes": reported → computed (subset)
    nexttile;
    n = numel(head_rep);
    idx = randperm(n, min(n, 80));        % at most 80 spokes for readability
    th_rep = deg2rad(head_rep(idx));
    th_est = deg2rad(ang_aligned(idx));
    polarplot([th_rep th_est].', [ones(numel(idx),1)*0.9 ones(numel(idx),1)*1.1].', ...
              '-', 'LineWidth', 0.7); hold on;
    polarplot(th_rep, ones(size(idx))*0.9, 'o', 'MarkerSize', 3, 'MarkerFaceColor',[.2 .2 .2], 'MarkerEdgeColor','none');
    polarplot(th_est, ones(size(idx))*1.1, 'o', 'MarkerSize', 3, 'MarkerFaceColor',[.1 .5 .9], 'MarkerEdgeColor','none');
    title('Reported → Computed (subset spokes)');

    % Print to console too
    fprintf('\n[plot_radflow_agreement]  MAE=%.2f°, med=%.2f° | ≤22.5°=%.1f%%, ≤11.25°=%.1f%% | ρ=%.3f\n', ...
            mae, med, pct22, pct11, rho);
end

% ---------- helpers ----------
function a = wrap360(a), a = mod(a,360); end
function a = wrap180(a), a = mod(a+180,360)-180; end

function [x,y] = wrapaware_xy(truth, est)
    % Make (x=truth, y=est) close to unity line by allowing y±360 shifts.
    T = wrap360(truth(:)); E = wrap360(est(:));
    % choose y, y+360, y-360 that is closest to x
    Ycands = [E, wrap360(E+360), wrap360(E-360)];
    [~, pick] = min(abs(Ycands - T), [], 2);
    y = Ycands(sub2ind(size(Ycands), (1:numel(E))', pick));
    x = T;
end

function ang_out = align_angles_to_truth(ang_in, truth_deg, is_centerout)
    % Try raw and +180° (for centerout), then subtract a single circular offset.
    C = [ang_in(:), wrap360(ang_in(:)+180)];
    if ~is_centerout, C = C(:,1); end
    best.mae = inf; best.ang = ang_in(:);
    for c = 1:size(C,2)
        a0 = C(:,c);
        delta = wrap180(a0 - truth_deg(:));
        off = atan2d(mean(sind(delta)), mean(cosd(delta)));  % circular mean offset
        a1  = wrap360(a0 - off);
        mae = mean(abs(wrap180(a1 - truth_deg(:))));
        if mae < best.mae, best.mae = mae; best.ang = a1; end
    end
    ang_out = best.ang;
end

function rho = circ_corr_deg(xdeg, ydeg)
    x = deg2rad(wrap360(xdeg(:))); y = deg2rad(wrap360(ydeg(:)));
    xm = atan2(mean(sin(x)),mean(cos(x)));
    ym = atan2(mean(sin(y)),mean(cos(y)));
    num = sum(sin(x - xm).*sin(y - ym));
    den = sqrt(sum(sin(x - xm).^2) * sum(sin(y - ym).^2));
    rho = num/den;
end
