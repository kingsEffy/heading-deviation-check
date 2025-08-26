function [ resA, resB, data]=validate_radflow_from_stc1_v2(matPath, stim_mode, condition)
% Validate your radial optic-flow pipeline on CRCNS stc-1 (MSTd/VIP).
% Runs two pipelines:
%   A) nearest-4: pick dataset columns nearest to [45,135,225,315]
%   B) fit-then-fold (RECOMMENDED): fit [Bx,By] on ALL headings, predict at 4 targets
%
% Usage:
%   validate_radflow_from_stc1_v2('MSTd.mat','sinkin','vis')
%   validate_radflow_from_stc1_v2('VIP.mat','centerout','vis')

if nargin < 2 || isempty(stim_mode),  stim_mode  = 'sinkin'; end
if nargin < 3 || isempty(condition),  condition  = 'vis';    end
assert(ismember(lower(stim_mode), {'sinkin','centerout'}), 'stim_mode must be sinkin|centerout');
assert(ismember(lower(condition), {'vis','ves'}), 'condition must be vis|ves');

%% 0) Load (using official field names)
S = load(matPath);
assert(isfield(S,'experiment1'),'Missing experiment1 in %s', matPath);
E = S.experiment1;
assert(isfield(E,'units') && ~isempty(E.units),'experiment1.units missing/empty');

cond = lower(condition);
U1 = E.units(1);
assert(isfield(U1,cond) && isfield(U1.(cond),'stim_global') && isfield(U1.(cond),'resp_global'), ...
    'Expected units(i).%s.stim_global/resp_global', cond);

theta_deg = wrap360(E.units(1).(cond).stim_global(:));  % K×1 headings (deg)
K = numel(theta_deg);
nUnits = numel(E.units);

Rfull = nan(nUnits,K);
for i = 1:nUnits
    ui = E.units(i).(cond);
    ti = wrap360(ui.stim_global(:));
    ri = ui.resp_global(:);
    assert(numel(ti)==K && numel(ri)==K,'Unit %d has mismatched heading/resp lengths', i);
    if any(abs(wrap180(ti - theta_deg))>1e-6), error('Unit %d headings differ from unit 1', i); end
    Rfull(i,:) = ri(:).';
end

fprintf('Loaded %s | cond=%s | nUnits=%d | nHeadings=%d\n', matPath, cond, nUnits, K);

% (Optional but recommended) Match your function’s normalization:
% linearfit_radical z-scores across directions per unit:contentReference[oaicite:1]{index=1}.
Rfull_z = zscore(Rfull, 0, 2);
data.theta_deg = theta_deg;   % 1×K headings (deg)
data.Rfull     = Rfull;       % nUnits×K trial-avg responses
data.Rfull_z   = Rfull_z;     % z-scored (what metrics used)
data.stim_mode = stim_mode;   % 'sinkin' or 'centerout'

dirlist = [45 135 225 315];
speeds  = [25 50 100 200]';                    % dummy; only slope column uses these
Sspd    = repmat(speeds', nUnits, 1);

%% A) nearest-4 (what you ran before)
idx4 = zeros(1,4);
for k = 1:4, [~,idx4(k)] = min(abs(wrap180(theta_deg - dirlist(k)))); end
R4_near = Rfull_z(:, idx4);

resA = run_linearfit_and_compare(R4_near, theta_deg, Rfull_z, dirlist, Sspd, speeds, stim_mode);
print_summary('Nearest-4', resA);

%% B) fit-then-fold (recommended)
th = deg2rad(theta_deg(:));            % K×1
X  = [cos(th) sin(th)];                % K×2 (global fit basis)
BxBy = zeros(nUnits,2);
for i = 1:nUnits
    y = Rfull_z(i,:).';
    B = X \ y;                         % 2×1 least-squares fit on ALL headings
    BxBy(i,:) = B.';
end
Xp = [cosd(dirlist(:)) sind(dirlist(:))];    % 4×2
R4_fit = (Xp * BxBy.').';                    % nUnits×4 predicted at 4 targets

resB = run_linearfit_and_compare(R4_fit, theta_deg, Rfull_z, dirlist, Sspd, speeds, stim_mode);
print_summary('Fit-then-fold', resB);

%% Side-by-side deltas
fprintf('\n=== Delta (Fit-then-fold minus Nearest-4) ===\n');
fprintf('ΔPop heading (deg)          : %+5.1f\n', angdiff(resB.pop_heading, resA.pop_heading));
fprintf('ΔCircular MAE (deg)         : %+5.2f (↓ is better)\n', resB.mae - resA.mae);
fprintf('ΔMedian error (deg)         : %+5.2f (↓ is better)\n', resB.med - resA.med);
fprintf('Δ%% within ±22.5° (pp)       : %+5.1f\n', resB.pct22 - resA.pct22);
fprintf('Δ%% within ±11.25° (pp)      : %+5.1f\n', resB.pct11 - resA.pct11);
fprintf('ΔCircular corr (rho)        : %+5.3f\n', resB.rho - resA.rho);

end % main

% ================= helpers =================
function out = run_linearfit_and_compare(R4, theta_deg, Rfull_z, dirlist, Sspd, speeds, stim_mode)
    % Run solver on R4 and compare to argmax headings.
    [~, ang_w, ~, quad_counts, ~, ~, UV] = ...
        linearfit_radical(R4, Sspd, dirlist, speeds, lower(stim_mode)); 
    ang_deg = wrap360(rad2deg(ang_w));

    % population heading from your UV summary
    U = UV(:,1); V = UV(:,2); U(~isfinite(U))=0; V(~isfinite(V))=0;
    pop_heading = wrap360(rad2deg(atan2(sum(V), sum(U))));

    % quadrant counts by angle bins (for reference)
    qb = angle_quadrants(ang_deg);

    % dataset-reported heading per unit = argmax over ALL headings
    [~, ixMax] = max(Rfull_z, [], 2);
    truth = wrap360(theta_deg(ixMax));

    % align your angles to truth (handle centerout 180° flip + single offset)
    ang_aligned = align_angles_to_truth(ang_deg, truth, strcmpi(stim_mode,'centerout'));

    err = abs(wrap180(ang_aligned - truth));
    out.mae  = mean(err);
    out.med  = median(err);
    out.pct22= mean(err<=22.5)*100;
    out.pct11= mean(err<=11.25)*100;
    out.rho  = circ_corr_deg(ang_aligned, truth);
    out.pop_heading = pop_heading;
    out.quad_sign   = quad_counts;
    out.quad_bins   = qb;

    out.ang_w   = ang_w;                                  % radians (per unit)
    out.ang_deg = mod(rad2deg(ang_w), 360);               % degrees (raw)
    out.UV      = UV;                                     % for quadrant plots if needed


end

function print_summary(tag, R)
    fprintf('\n=== %s (%s geometry) ===\n', tag, 'linearfit');
    fprintf('Population heading (deg): %5.1f\n', R.pop_heading);
    fprintf('Quadrants by sign  [RU;LU;LD;RD] = [%d;%d;%d;%d]\n', R.quad_sign);
    fprintf('Quadrants by angle [RU,LU,LD,RD] = [%d,%d,%d,%d]\n', R.quad_bins);
    fprintf('Circular MAE = %5.2f°, median = %5.2f° |  within ±22.5°: %4.1f%%, ±11.25°: %4.1f%% | rho = %.3f\n', ...
        R.mae, R.med, R.pct22, R.pct11, R.rho);
end

function qb = angle_quadrants(ang_deg)
    a = wrap360(ang_deg);
    qb = [ sum(a>=315 | a<45); sum(a>=45  & a<135); sum(a>=135 & a<225); sum(a>=225 & a<315) ];
end

function a = wrap360(a), a = mod(a,360); end
function a = wrap180(a), a = mod(a+180,360)-180; end
function d = angdiff(a,b), d = wrap180(a-b); end

function ang_out = align_angles_to_truth(ang_in, truth_deg, is_centerout)
    C = [ang_in(:), wrap360(ang_in(:)+180)];
    if ~is_centerout, C = C(:,1); end
    best.mae = inf; best.ang = ang_in(:);
    for c = 1:size(C,2)
        a0 = C(:,c);
        delta = wrap180(a0 - truth_deg(:));
        off = atan2d(mean(sind(delta)), mean(cosd(delta)));  % circular mean offset
        a1  = wrap360(a0 - off);
        err = abs(wrap180(a1 - truth_deg(:)));
        m   = mean(err);
        if m < best.mae, best.mae = m; best.ang = a1; end
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

