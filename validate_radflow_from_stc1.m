function validate_radflow_from_stc1(matPath, stim_mode, condition)
% Validate radial optic-flow pipeline on CRCNS stc-1 (MSTd/VIP) data.
% 
% validate_radflow_from_stc1('C:\zy\code\mcode\adofflineanalysis-main\MSTd.mat', 'sinkin', 'vis')

% Usage:
%   validate_radflow_from_stc1('VIP.mat', 'sinkin', 'vis')
%   validate_radflow_from_stc1('MSTd.mat','centerout','vis')
%
% - matPath    : path to VIP.mat or MSTd.mat from CRCNS stc-1
% - stim_mode  : 'sinkin' or 'centerout' (your geometry convention)
% - condition  : 'vis' or 'ves' (visual vs vestibular condition in dataset)
%
% This follows the official stc-1 demo field names:
%   experiment1.units(i).vis.stim_global, .resp_global   (and .ves.*)  [VIP.mat]
% and keeps your linearfit_radical UNCHANGED.
%
% Reference for loading fields: VIP.mat demo code (CRCNS stc-1)  << IMPORTANT
% (load VIP.mat; units=experiment1.units; plot units(i).vis.stim_global vs resp_global)
% Source: the demo you attached.  (We mirror that field structure.) 
%
% Outputs printed to console:
%   - Population heading (deg) from your UV (suminf-safe)
%   - Quadrant counts [RU; LU; LD; RD] from your sign rules
%   - Circular error metrics vs dataset-reported per-unit headings (argmax)
%   - Agreement between sign-based quadrants and angle-binned quadrants

    if nargin < 2 || isempty(stim_mode),  stim_mode  = 'sinkin';   end
    if nargin < 3 || isempty(condition),  condition  = 'vis';       end
    assert(ismember(lower(stim_mode), {'sinkin','centerout'}), 'stim_mode must be sinkin|centerout');
    assert(ismember(lower(condition), {'vis','ves'}), 'condition must be vis|ves');

    %% ---------- 0) Load dataset (VIP.mat or MSTd.mat) ----------
    S = load(matPath);
    assert(isfield(S,'experiment1'), 'Expected field experiment1 is missing in %s', matPath);
    E = S.experiment1;  % global headings block

    % Official stc-1 demo uses fields: units(i).vis.stim_global / resp_global (or .ves.*)
    % We'll adhere to that (no guessing) based on your attached example.
    % If the file has no 'vis' (or 'ves'), this will error clearly.
    cond = lower(condition);
    assert(isfield(E,'units') && ~isempty(E.units), 'experiment1.units missing or empty.');
    U1 = E.units(1);
    assert(isfield(U1, cond), 'units(1).%s not found in %s', cond, matPath);
    assert(isfield(U1.(cond),'stim_global') && isfield(U1.(cond),'resp_global'), ...
        'units(1).%s.stim_global/resp_global not found (check file).', cond);

    % Heading vector (deg), consistent across units
    theta_deg = wrap360(E.units(1).(cond).stim_global(:));
    K         = numel(theta_deg);
    assert(K >= 4, 'Need at least 4 global headings.');

    % Build Rfull: [nUnits x K] trial-averaged responses for the chosen condition
    nUnits = numel(E.units);
    Rfull  = nan(nUnits, K);
    for i = 1:nUnits
        ui = E.units(i).(cond);
        ti = wrap360(ui.stim_global(:));
        ri = ui.resp_global(:);
        assert(numel(ti)==K && numel(ri)==K, ...
            'Unit %d has inconsistent stim/resp length (%d vs %d).', i, numel(ti), numel(ri));
        % Ensure alignment (should already match)
        if any(abs(wrap180(ti - theta_deg)) > 1e-6)
            error('Unit %d has different heading vector than units(1).', i);
        end
        Rfull(i,:) = ri(:).';
    end

    fprintf('Loaded %s | condition=%s | nUnits=%d | nHeadings=%d\n', matPath, cond, nUnits, K);

    %% ---------- 1) Fold to your 4 radial headings [45 135 225 315] ----------
    % dirlist = [45 135 225 315];                 % your basis for linearfit_radical
    % idx4 = zeros(1,4);
    % for k = 1:4
    %     [~, idx4(k)] = min(abs(wrap180(theta_deg - dirlist(k))));
    % end
    % R4 = Rfull(:, idx4);                        % nUnits x 4 (trial-avg visual responses at nearest radial)


% --- Replace "nearest 4 columns" with cosine fit over ALL headings ---
th = deg2rad(theta_deg(:));                 % K×1
X  = [cos(th) sin(th)];                     % K×2

BxBy = zeros(size(Rfull,1),2);              % nUnits×2
for i = 1:size(Rfull,1)
    y = Rfull(i,:).';                       % K×1 (trial-avg visual resp)
    % (Optional) zscore across headings to match your function’s internal normalization:
    y = (y - mean(y)) / std(y);             % improves robustness
    B  = X \ y;                              % 2×1 least-squares
    BxBy(i,:) = B.';
end

% Predict each unit’s responses at your 4 radial targets
dirlist = [45 135 225 315];
Xp = [cosd(dirlist(:)) sind(dirlist(:))];   % 4×2
R4  = (Xp * BxBy.').';                      % nUnits×4   (predicted responses at 4 targets)



    %% ---------- 2) Run your solver (UNCHANGED) ----------
    % speeds only fill the slope column; they don’t affect angle/quadrant results.
    speeds = [25 50 100 200]';                  % dummy (deg/s)
    Sspd   = repmat(speeds', nUnits, 1);

    [flow_strength, ang_w, T, quad_counts, labels, P, UV] = ...
        linearfit_radical(R4, Sspd, dirlist, speeds, lower(stim_mode)); %#ok<ASGLU>

    % Population heading from your UV (your convention)
    U = UV(:,1); V = UV(:,2);
    U(~isfinite(U)) = 0;  V(~isfinite(V)) = 0;
    pop_heading = wrap360(rad2deg(atan2(sum(V), sum(U))));

    % Per-unit angles from your solver
    ang_deg = wrap360(rad2deg(ang_w));

    % Angle-binned quadrant counts
    [ruA, luA, ldA, rdA] = bucket4(ang_deg);

    fprintf('\n=== Linearfit summary (%s) ===\n', lower(stim_mode));
    fprintf('Population heading (deg): %.1f\n', pop_heading);
    fprintf('Quadrants by sign [RU;LU;LD;RD] = [%d;%d;%d;%d]\n', quad_counts);
    fprintf('Quadrants by angle bins [RU,LU,LD,RD] = [%d,%d,%d,%d]\n', ruA,luA,ldA,rdA);

    %% ---------- 3) Confirm vs the dataset-reported “heading” per unit ----------
    % Reported heading = argmax across ALL global headings (richer than 4-col fold)
    [~, ixMax]  = max(Rfull, [], 2);
    head_rep    = wrap360(theta_deg(ixMax));

    % Align your per-unit angles to reported headings (centerout flip + single offset)
    ang_aligned = align_angles_to_truth(ang_deg, head_rep, strcmpi(stim_mode,'centerout'));

    % Circular error metrics
    err = abs(wrap180(ang_aligned - head_rep));
    mae = mean(err); med = median(err);
    pct22 = mean(err <= 22.5)*100;
    pct11 = mean(err <= 11.25)*100;
    rho = circ_corr_deg(ang_aligned, head_rep);

    fprintf('\n=== Heading agreement (your angle vs argmax heading) ===\n');
    fprintf('Circular MAE = %.2f°, median = %.2f°\n', mae, med);
    fprintf('%% within ±22.5°: %.1f%%   |  %% within ±11.25°: %.1f%%\n', pct22, pct11);
    fprintf('Circular corr (Fisher–Lee rho): %.3f\n', rho);

    % Quadrant agreement (your sign vs aligned angles)
    [ruB, luB, ldB, rdB] = bucket4(ang_aligned);
    fprintf('Quadrants by aligned angle bins [RU,LU,LD,RD] = [%d,%d,%d,%d]\n', ruB,luB,ldB,rdB);
end

%% ======================= helpers (local) =========================
function a = wrap360(a), a = mod(a,360); end
function a = wrap180(a), a = mod(a+180,360)-180; end

function [ru,lu,ld,rd] = bucket4(ang_deg)
    a = wrap360(ang_deg);
    ru = sum(a>=315 | a<45);
    lu = sum(a>=45  & a<135);
    ld = sum(a>=135 & a<225);
    rd = sum(a>=225 & a<315);
end

function ang_out = align_angles_to_truth(ang_in, truth_deg, is_centerout)
    % Try your raw angles and +180° (for centerout basis), then remove a single
    % global circular offset that minimises median circular error.
    C = [ang_in(:), wrap360(ang_in(:)+180)];
    if ~is_centerout, C = C(:,1); end
    best.mae = inf; best.ang = ang_in(:);

    for c = 1:size(C,2)
        a0 = C(:,c);
        delta = wrap180(a0 - truth_deg(:));
        % circular mean offset (deg)
        off = atan2d(mean(sind(delta)), mean(cosd(delta)));
        a1  = wrap360(a0 - off);
        err = abs(wrap180(a1 - truth_deg(:)));
        mae = mean(err);
        if mae < best.mae
            best.mae = mae; best.ang = a1;
        end
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
