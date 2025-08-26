function [res, data] = validate_RSI_from_stc1(matPath, stim_mode, condition, do_plot)
% Validate RSI_radicalflow on CRCNS stc-1 (MSTd/VIP).
%
% Usage:
%   [res, data] = validate_RSI_from_stc1('MSTd.mat','sinkin','vis', true);
%   [res, data] = validate_RSI_from_stc1('VIP.mat','centerout','vis', true);
%
% Inputs
%   matPath   : VIP.mat or MSTd.mat (CRCNS stc-1)
%   stim_mode : 'sinkin' or 'centerout' (passed to RSI_radicalflow)
%   condition : 'vis' or 'ves'
%   do_plot   : true/false – if true, calls plot_radflow_agreement(...)
%
% Outputs
%   res.best_deg      : per-unit RSI preferred angle (deg, 0:45:315)
%   res.pop_heading   : RSI population heading (deg)
%   res.metrics       : struct with MAE/median/%≤22.5/%≤11.25/ρ
%   data.theta_deg    : dataset headings (deg) used for argmax truth
%   data.Rfull_z      : z-scored nUnits×K responses used for truth
%   data.R8           : nUnits×8 matrix fed to RSI_radicalflow
%
% Notes
%   - We fit [Bx,By] on ALL dataset headings (K=10 in MSTd), then predict at
%     RSI's 8 directions (0:45:315) to avoid quantization/misalignment.
%   - We keep your RSI function untouched and compare RSI's per-unit angles
%     to “reported” headings (argmax over ALL headings, z-scored).
%
% Requires: your RSI_radicalflow.m, and (optional) plot_radflow_agreement.m

if nargin < 2 || isempty(stim_mode),  stim_mode  = 'sinkin'; end
if nargin < 3 || isempty(condition),  condition  = 'vis';    end
if nargin < 4, do_plot = false; end
assert(ismember(lower(stim_mode),{'sinkin','centerout'}), 'stim_mode must be sinkin|centerout');
assert(ismember(lower(condition),{'vis','ves'}), 'condition must be vis|ves');

%% 0) Load dataset (official field names)
S = load(matPath);
assert(isfield(S,'experiment1'), 'experiment1 missing from %s', matPath);
E = S.experiment1;
assert(isfield(E,'units') && ~isempty(E.units), 'experiment1.units missing/empty');

cond = lower(condition);
U1 = E.units(1);
assert(isfield(U1,cond) && isfield(U1.(cond),'stim_global') && isfield(U1.(cond),'resp_global'), ...
    'Expected units(i).%s.stim_global/resp_global fields not found.', cond);

theta_deg = wrap360(E.units(1).(cond).stim_global(:));     % K×1 headings (deg)
K = numel(theta_deg);
nUnits = numel(E.units);

Rfull = nan(nUnits, K);
for i = 1:nUnits
    ui = E.units(i).(cond);
    ti = wrap360(ui.stim_global(:));
    ri = ui.resp_global(:);
    assert(numel(ti)==K && numel(ri)==K, 'Unit %d heading/resp length mismatch', i);
    if any(abs(wrap180(ti - theta_deg))>1e-6), error('Unit %d has different heading vector', i); end
    Rfull(i,:) = ri(:).';
end

fprintf('Loaded %s | cond=%s | nUnits=%d | nHeadings=%d\n', matPath, cond, nUnits, K);

% Match your previous validation: z-score across headings per unit (shape)
Rfull_z = zscore(Rfull, 0, 2);

%% 1) Build the 8-direction matrix for RSI by fitting ALL headings
angles8 = (0:45:315).';                   % 8×1
Xp8     = [cosd(angles8) sind(angles8)];  % 8×2
th      = deg2rad(theta_deg(:));          % K×1
X       = [cos(th) sin(th)];              % K×2

BxBy = zeros(nUnits,2);
for i = 1:nUnits
    y = Rfull_z(i,:).';
    B = X \ y;                            % 2×1 least-squares on ALL K headings
    BxBy(i,:) = B.';
end
R8 = (Xp8 * BxBy.').';                    % nUnits×8 predicted at RSI angles

%% 2) Run RSI_radicalflow (UNCHANGED)
resolution = [4320,4320];% [1024 1024];                 % square to avoid ellipse bias
[RSI_ik, best_deg, RSI_global, pop_vec, fwd, dev] = ...
    RSI_radicalflow(R8, resolution, lower(stim_mode));      %#ok<ASGLU>

% population heading from RSI:
pop_heading = wrap360(rad2deg(atan2(pop_vec(2), pop_vec(1))));

%% 3) “Reported” heading per unit = argmax over ALL headings (truth)
[~, ixMax] = max(Rfull_z, [], 2);
truth = wrap360(theta_deg(ixMax));

% Align RSI angles to truth (handle centerout 180° & single global offset)
best_aligned = align_angles_to_truth(best_deg, truth, strcmpi(stim_mode,'centerout'));

% Circular errors
err   = wrap180(best_aligned - truth);
absE  = abs(err);

metrics.MAE   = mean(absE);
metrics.med   = median(absE);
metrics.pct22 = mean(absE <= 22.5)*100;
metrics.pct11 = mean(absE <= 11.25)*100;
metrics.rho   = circ_corr_deg(best_aligned, truth);

%% 4) Package results
res.best_deg    = best_deg(:);
res.pop_heading = pop_heading;
res.metrics     = metrics;

data.theta_deg  = theta_deg(:).';
data.Rfull_z    = Rfull_z;
data.R8         = R8;
data.stim_mode  = stim_mode;

fprintf('\n=== RSI validation (%s / %s) ===\n', upper(condition), lower(stim_mode));
fprintf('Population heading (deg): %.1f\n', res.pop_heading);
fprintf('MAE = %.2f°, median = %.2f° |  ≤22.5° = %.1f%%, ≤11.25° = %.1f%% | ρ = %.3f\n', ...
        metrics.MAE, metrics.med, metrics.pct22, metrics.pct11, metrics.rho);

%% 5) Optional plots (reuse your linearfit figure helper)
if do_plot
    try
        plot_radflow_agreement(theta_deg, Rfull_z, deg2rad(best_deg), lower(stim_mode), ...
            sprintf('RSI — %s (%s)', upper(condition), lower(stim_mode)));
    catch ME
        warning('plot_radflow_agreement failed: %s\n(Ensure it is on the path.)', ME.message);
    end
end
end

% ================= helpers =================
function a = wrap360(a), a = mod(a,360); end
function a = wrap180(a), a = mod(a+180,360)-180; end

function ang_out = align_angles_to_truth(ang_in_deg, truth_deg, is_centerout)
    % Try raw and +180° (for centerout), then subtract a single circular offset.
    C = [ang_in_deg(:), wrap360(ang_in_deg(:)+180)];
    if ~is_centerout, C = C(:,1); end
    best.mae = inf; best.ang = ang_in_deg(:);
    for c = 1:size(C,2)
        a0 = C(:,c);
        delta = wrap180(a0 - truth_deg(:));
        off = atan2d(mean(sind(delta)), mean(cosd(delta)));   % circular mean offset
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
