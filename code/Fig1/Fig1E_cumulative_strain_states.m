%% SingleExperiment_X_CumulativeStrainStates_dprime15_Frame0_Frame1.m
% ONE experiment ONLY (X condition):
% - ROI parsed from experiment .txt "light box" (pixel coords)
% - Load ALL PIV frames in PIV_Results
% - Compute strain field per frame from u(x,y), v(x,y) (um per frame)
% - Build CUMULATIVE strain field: C_k(x,y) = sum_{j<=k} F_j(x,y)
% - For each frame k: extract ROI distribution of C_k values on the PIV grid
%
% STATES (your request):
%   State 1 = frame 0 baseline: C_0 = 0 everywhere  (conceptual, no PIV file)
%   State 2 = frame 1 (the first processed PIV frame in this script)
%   State 3+:
%       Greedy rule: pick earliest later frame whose ROI distribution has
%       |d'| >= DPRIME_THRESH relative to previous picked state's ROI distribution.
%
% Notes:
% - "frame 1" here means the first entry in the PIV list after SKIP_FIRST_N and MAX_FRAMES.
% - d' is computed from spatial samples inside ROI (correlated); operational criterion.

clearvars; close all; clc;

%% ================== USER SETTINGS ==================
SET_ID = 3;   % <-- choose one X experiment

BASE_DIR_X = 'PATH_TO_DATA/Experiments/SingleLightSpot_Data&PIV/';

% ---- PIV file convention ----
PIV_STEP      = 2;      % PIV_ti - t(i+PIV_STEP)
POS_COL_X_PX  = 1;      % x position (pixels) in PIV txt
POS_COL_Y_PX  = 2;      % y position (pixels)
U_COL         = 13;     % ux0 column
V_COL         = 14;     % uy0 column

% Units (used if parsing fails)
DT_SEC = 1/3;   % seconds per PIV frame
DEFAULT_X_UNIT = 0.1;   % um/px
DEFAULT_Y_UNIT = 0.1;   % um/px

% ---- per-frame strain field F to accumulate ----
% 'epsVol'    : divergence = exx + eyy
% 'epsVolAbs' : abs(divergence)
% 'epsDev'    : deviatoric strain magnitude
FIELD_METRIC = 'epsVolAbs';

% ---- strain computation knobs ----
SMOOTH_UV_SIGMA_GRID = 0.0;   % in GRID units; try 0.8–1.2 if derivatives are noisy
USE_MEDIAN_FILL      = true;  % fill NaNs by local median before smoothing/gradients

% ---- time/frame handling ----
SKIP_FIRST_N = 0;            % keep 0 to make "frame1" be the earliest PIV frame
MAX_FRAMES   = inf;          % cap frames processed (inf for all)

% ---- d' threshold for state separation (state3+) ----
DPRIME_THRESH = 2;

% Guards
MIN_ROI_PTS   = 50;          % require at least this many ROI samples
REQUIRE_MONOTONIC_MEAN = true; % keep states increasing in mean cumulative metric
ALLOW_EQUAL_MEAN = true;

% ---- output ----
SAVE_DIR = fullfile(pwd, 'SingleExp_X_CumulativeStrainStates');
SAVE_FIG = true;
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end

%% ================== PATHS ==================
expDir = fullfile(BASE_DIR_X, sprintf('set%d', SET_ID));
pivDir = fullfile(expDir, 'PIV_Results');

if ~isfolder(expDir), error('Missing expDir: %s', expDir); end
if ~isfolder(pivDir), error('Missing pivDir: %s', pivDir); end

%% ================== PARSE UNITS ==================
[x_unit, y_unit] = parse_xy_units_from_any_txt(expDir);
if isnan(x_unit), x_unit = DEFAULT_X_UNIT; end
if isnan(y_unit), y_unit = DEFAULT_Y_UNIT; end

%% ================== ROI (X only: from txt pixels) ==================
[center_px, radius_px] = parse_roi_from_txt_pixels(expDir, false, SET_ID);
if any(isnan([center_px, radius_px]))
    error('Failed ROI parse from txt for set%d.', SET_ID);
end
center_um = [center_px(1) * x_unit, center_px(2) * y_unit];
radius_um = radius_px * sqrt(x_unit * y_unit); % effective radius in um

fprintf('set%d: x_unit=%.4g um/px, y_unit=%.4g um/px\n', SET_ID, x_unit, y_unit);
fprintf('ROI center(px)=(%.2f,%.2f), radius(px)=%.2f; center(um)=(%.2f,%.2f), radius(um)=%.2f\n', ...
    center_px(1), center_px(2), radius_px, center_um(1), center_um(2), radius_um);

%% ================== LIST PIV FRAMES ==================
i_list = list_all_piv_i(pivDir, PIV_STEP);
if isempty(i_list), error('No PIV files found in %s', pivDir); end
i_list = i_list(:)';

if SKIP_FIRST_N > 0
    i_list = i_list(min(SKIP_FIRST_N+1, numel(i_list)):end);
end
if isfinite(MAX_FRAMES)
    i_list = i_list(1:min(numel(i_list), MAX_FRAMES));
end

nF = numel(i_list);
fprintf('Processing %d PIV frames\n', nF);
if nF < 1, error('No frames after trimming.'); end

%% ================== COMPUTE CUMULATIVE FIELD + ROI DISTRIBUTIONS ==================

tIdx = (1:nF)';                 % frame index (keep this if needed)
tSec = (tIdx - 1) * DT_SEC;     % time in seconds

cumMean = nan(nF,1);
cumStd  = nan(nF,1);
roiN    = nan(nF,1);
frameOK = false(nF,1);
roiValsPerFrame = cell(nF,1);

C = [];
roiMask = [];
gridSet = false;

for k = 1:nF
    i = i_list(k);

    [Xum, Yum, Ug, Vg, ok] = load_uv_grid_um(pivDir, i, PIV_STEP, ...
        POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, x_unit, y_unit);
    if ~ok, continue; end

    [epsVol, epsDev, okS] = compute_strain_from_uv(Xum, Yum, Ug, Vg, ...
        USE_MEDIAN_FILL, SMOOTH_UV_SIGMA_GRID);
    if ~okS, continue; end

    switch lower(FIELD_METRIC)
        case 'epsvol'
            F = epsVol;
        case 'epsvolabs'
            F = abs(epsVol);
        case 'epsdev'
            F = epsDev;
        otherwise
            error('Unknown FIELD_METRIC: %s', FIELD_METRIC);
    end

    if ~gridSet
        C = zeros(size(F), 'like', F);
        roiMask = hypot(Xum - center_um(1), Yum - center_um(2)) <= radius_um;
        gridSet = true;
    else
        if ~isequal(size(F), size(C))
            warning('Frame k=%d has different grid size; skipping.', k);
            continue;
        end
    end

    C = C + F;

    vals = C(roiMask);
    vals = vals(isfinite(vals));

    roiN(k) = numel(vals);
    if roiN(k) < MIN_ROI_PTS
        continue;
    end

    roiValsPerFrame{k} = vals;
    cumMean(k) = mean(vals);
    cumStd(k)  = std(vals, 0);
    frameOK(k) = true;
end

validIdx = find(frameOK);
if isempty(validIdx)
    error('No valid frames (MIN_ROI_PTS=%d).', MIN_ROI_PTS);
end

%% ================== DEFINE STATE1 AND STATE2 ==================
% State1: frame0 baseline (conceptual, not in data arrays)
% State2: MUST be frame1 (k=1) in this script.
state2_k = 1;

if ~frameOK(state2_k)
    % If frame 1 is invalid, we still obey your rule conceptually, but we
    % must have a distribution to compute d' for later states. So we hard-fail.
    error('State2 must be frame1 (k=1), but frame1 is invalid (not enough ROI points or failed load). Fix data/parameters.');
end

picked = [];
picked(end+1) = state2_k; %#ok<AGROW>
prevK = state2_k;
prevVals = roiValsPerFrame{prevK};
prevMu = cumMean(prevK);

%% ================== GREEDY PICK STATE3+ USING d' THRESH ==================
for k = prevK+1:nF
    if ~frameOK(k), continue; end

    currVals = roiValsPerFrame{k};
    if isempty(currVals) || numel(currVals) < MIN_ROI_PTS
        continue;
    end

    if REQUIRE_MONOTONIC_MEAN
        mu = cumMean(k);
        if ~isfinite(mu), continue; end
        if ALLOW_EQUAL_MEAN
            if mu < prevMu, continue; end
        else
            if mu <= prevMu, continue; end
        end
    end

    dp = dprime_from_samples(currVals, prevVals);

    if isfinite(dp) && abs(dp) >= DPRIME_THRESH
        picked(end+1) = k; %#ok<AGROW>
        prevK = k;
        prevVals = currVals;
        prevMu = cumMean(k);
    end
end

totalStates = numel(picked) + 1; % +1 for baseline state1
bits = log2(max(totalStates,1));

fprintf('\n===== STATE PICKING (state1=frame0, state2=frame1, d''>=%.2f for state3+) =====\n', DPRIME_THRESH);
fprintf('State1: baseline frame0 (C=0)\n');
fprintf('State2: frame1 in script (k=1) -> PIV i=%d\n', i_list(1));
fprintf('Metric accumulated: %s\n', FIELD_METRIC);
fprintf('Picked %d additional states after state2 => total states=%d, bits=%.3f\n', ...
    numel(picked)-1, totalStates, bits);
fprintf('Picked k indices [state2..]: %s\n', mat2str(picked(:)'));
fprintf('Picked PIV i-values: %s\n', mat2str(i_list(picked)));

%% ================== COMPUTE d' BETWEEN PICKED STATES (2-3, 3-4, ...) ==================
if numel(picked) < 2
    fprintf('No adjacent picked-state pairs beyond state2.\n');
    dprime_pairs = [];
else
    dprime_pairs = nan(numel(picked)-1,1);
    for s = 1:(numel(picked)-1)
        kA = picked(s);
        kB = picked(s+1);
        A = roiValsPerFrame{kA};
        B = roiValsPerFrame{kB};
        dprime_pairs(s) = dprime_from_samples(B, A); % next vs prev
        fprintf('d''(state%d vs state%d) = %.3f  (k=%d vs %d)\n', ...
            s+1, s+2, dprime_pairs(s), kA, kB);
    end
    fprintf('Min |d''| among adjacent picked states (2+) = %.3f\n', min(abs(dprime_pairs), [], 'omitnan'));
end

%% ================== PLOT 1: ROI mean ± ROI std and picked states ==================
% ---- color controls ----
grayLevel   = 0.6;                     % 0 = black, 1 = white (tune this)
rgb_Gray    = grayLevel * [1 1 1];

rgb_Black   = [0 0 0];
rgb_Purple  = [135, 85, 155] / 255;


% ---- styling knobs ----
FZ = 20;
FZ_LEG = 15;
LW = 2;
MS = 4;
MS_PICK = 9;

tSec = (tIdx - 1) * DT_SEC;

fig1 = figure('Color','w'); hold on;

valid = frameOK & isfinite(cumMean) & isfinite(cumStd);

% Main errorbar curve (gray)
eh = errorbar(tSec(valid), cumMean(valid), cumStd(valid), 'o-', ...
    'LineWidth', LW, 'MarkerSize', MS, 'Color', rgb_Gray);
eh.CapSize = 6;

% Baseline (state1 = 0) in black
hBase = yline(0, '--', 'LineWidth', LW, 'Color', rgb_Black);

% Dummy handle for legend (picked states)
hPickLegend = plot(nan, nan, 'o', ...
    'MarkerSize', MS_PICK, 'LineWidth', LW, 'Color', rgb_Purple);

% Picked states + labels (purple)
for s = 1:numel(picked)
    k = picked(s);

    plot(tSec(k), cumMean(k), 'o', ...
        'MarkerSize', MS_PICK, ...
        'LineWidth', LW, ...
        'Color', rgb_Purple);

    text(tSec(k), cumMean(k), sprintf('  state%d', s+1), ...
        'FontSize', 12, ...
        'VerticalAlignment','bottom', ...
        'Color', rgb_Purple);
end

xlabel(sprintf('Time (s)', DT_SEC), 'FontSize', FZ);
ylabel('Cumulative strain', 'FontSize', FZ);

leg1 = legend([eh, hBase, hPickLegend], ...
    {'ROI mean \pm ROI std', 'Baseline state1 (=0)', 'Picked states'}, ...
    'Location','southeast');
set(leg1, 'Box','on', 'FontSize', FZ_LEG-3);

set(gca, 'FontSize', FZ, 'LineWidth', 1, 'TickLength', [0.02 0.02]);
grid on; box on;

if SAVE_FIG
    outP = fullfile(SAVE_DIR, sprintf( ...
        'set%d_X_cumulative_%s_states_dprime%.1f_time.eps', ...
        SET_ID, FIELD_METRIC, DPRIME_THRESH));

    exportgraphics(fig1, outP, ...
        'ContentType', 'vector', ...   % <-- critical for EPS
        'Resolution', 300);            % resolution ignored for pure vector but safe

    fprintf('Saved figure: %s\n', outP);
end


%% ================== SAVE TABLE ==================
outT = table();
outT.FrameIdx = tIdx;
outT.PIV_i    = i_list(:);
outT.OK       = frameOK(:);
outT.ROI_N    = roiN(:);
outT.CumMean  = cumMean(:);
outT.CumStd   = cumStd(:);

outT.Picked   = false(nF,1);
outT.Picked(picked) = true;

% State numbering: 0=none, 2.. for picked states (state1 is conceptual baseline)
stateNum = zeros(nF,1);
for s = 1:numel(picked)
    stateNum(picked(s)) = s+1; % state2=2, state3=3,...
end
outT.StateNum = stateNum;

outCSV = fullfile(SAVE_DIR, sprintf('set%d_X_cumulative_%s_timeseries_dprime%.1f_f0f1.csv', SET_ID, FIELD_METRIC, DPRIME_THRESH));
writetable(outT, outCSV);
fprintf('Saved timeseries CSV: %s\n', outCSV);

fprintf('Done.\n');

%% ================== FUNCTIONS ==================

function dp = dprime_from_samples(A, B)
    A = A(isfinite(A));
    B = B(isfinite(B));
    if numel(A) < 2 || numel(B) < 2
        dp = nan; return;
    end
    muA = mean(A); muB = mean(B);
    sA  = std(A, 0); sB = std(B, 0);
    denom = sqrt(0.5*(sA^2 + sB^2));
    dp = (muA - muB) / max(denom, eps);
end

function [Xum, Yum, Ug, Vg, ok] = load_uv_grid_um(pivDir, i, pivStep, ...
                                                 POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
                                                 x_unit, y_unit)
    ok = false;
    Xum = []; Yum = []; Ug = []; Vg = [];

    pivFile = fullfile(pivDir, sprintf('PIV_t%d-t%d.txt', i, i+pivStep));
    if exist(pivFile, 'file') ~= 2, return; end

    data = load(pivFile);
    needCols = max([POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL]);
    if size(data,2) < needCols, return; end

    x_um = x_unit * data(:, POS_COL_X_PX);
    y_um = y_unit * data(:, POS_COL_Y_PX);

    u_um = x_unit * data(:, U_COL) / pivStep;
    v_um = y_unit * data(:, V_COL) / pivStep;

    [Xum, Yum, Ug, Vg] = vector_to_grid(x_um, y_um, u_um, v_um);
    ok = ~isempty(Xum) && all(size(Xum) == size(Ug));
end

function [epsVol, epsDev, ok] = compute_strain_from_uv(Xum, Yum, Ug, Vg, ...
                                                       doMedianFill, smoothSigma)
    ok = false;
    epsVol = []; epsDev = [];

    if isempty(Ug) || isempty(Vg), return; end

    U = Ug; V = Vg;

    if doMedianFill
        U = fill_nans_local(U, 2);
        V = fill_nans_local(V, 2);
    end

    if smoothSigma > 0
        U = imgaussfilt(U, smoothSigma, 'Padding','replicate');
        V = imgaussfilt(V, smoothSigma, 'Padding','replicate');
    end

    dx = median(diff(Xum(1,:)));
    dy = median(diff(Yum(:,1)));
    if ~isfinite(dx) || ~isfinite(dy) || dx <= 0 || dy <= 0, return; end

    [dUdy, dUdx] = gradient(U, dy, dx);
    [dVdy, dVdx] = gradient(V, dy, dx);

    exx = dUdx;
    eyy = dVdy;
    exy = 0.5*(dUdy + dVdx);

    epsVol = exx + eyy;
    epsDev = sqrt((exx - eyy).^2 + (2*exy).^2) / 2;

    ok = true;
end

function A = fill_nans_local(A, rad)
    nanMask = ~isfinite(A);
    if ~any(nanMask(:)), return; end

    [ny, nx] = size(A);
    B = A;

    for iy = 1:ny
        y1 = max(1, iy-rad); y2 = min(ny, iy+rad);
        for ix = 1:nx
            if ~nanMask(iy,ix), continue; end
            x1 = max(1, ix-rad); x2 = min(nx, ix+rad);
            win = A(y1:y2, x1:x2);
            win = win(isfinite(win));
            if isempty(win)
                B(iy,ix) = 0;
            else
                B(iy,ix) = median(win);
            end
        end
    end
    A = B;
end

function [x_unit, y_unit] = parse_xy_units_from_any_txt(expDir)
    x_unit = nan; y_unit = nan;
    txtFiles = dir(fullfile(expDir, '*.txt'));
    txtFiles = txtFiles(~startsWith({txtFiles.name}, '._'));
    if isempty(txtFiles), return; end

    for k = 1:numel(txtFiles)
        fp = fullfile(expDir, txtFiles(k).name);
        try
            lines = readlines(fp);
        catch
            continue;
        end
        if numel(lines) < 6, continue; end

        xu = parse_last_number(lines(5));
        yu = parse_last_number(lines(6));
        if ~isnan(xu), x_unit = xu; end
        if ~isnan(yu), y_unit = yu; end
        if ~isnan(x_unit) && ~isnan(y_unit), return; end
    end
end

function val = parse_last_number(lineStr)
    s = char(lineStr);
    nums = regexp(s, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
    if isempty(nums), val = nan; else, val = str2double(nums{end}); end
end

function i_list = list_all_piv_i(pivDir, pivStep)
    files = dir(fullfile(pivDir, 'PIV_t*-t*.txt'));
    ii = [];
    for k = 1:numel(files)
        tok = regexp(files(k).name, 'PIV_t(\d+)-t(\d+)\.txt', 'tokens', 'once');
        if isempty(tok), continue; end
        i1 = str2double(tok{1});
        i2 = str2double(tok{2});
        if ~isnan(i1) && ~isnan(i2) && (i2 == i1 + pivStep)
            ii(end+1) = i1; %#ok<AGROW>
        end
    end
    i_list = sort(unique(ii));
end

function [center_px, radius_px] = parse_roi_from_txt_pixels(expDir, dbg, setID)
    center_px = [nan nan];
    radius_px = nan;

    txtFiles = dir(fullfile(expDir, '*.txt'));
    txtFiles = txtFiles(~startsWith({txtFiles.name}, '._'));
    if isempty(txtFiles)
        if dbg, fprintf('set%d: no txt files.\n', setID); end
        return;
    end

    for k = 1:numel(txtFiles)
        fp = fullfile(expDir, txtFiles(k).name);
        try
            lines = readlines(fp);
        catch
            continue;
        end

        [bestIdx, pts] = find_best_point_line(lines);
        if dbg && ~isnan(bestIdx)
            fprintf('set%d ROI line %d: %s\n', setID, bestIdx, char(lines(bestIdx)));
        end

        if size(pts,1) < 2, continue; end

        p1 = pts(1,:); p2 = pts(2,:);
        center_px = 0.5*(p1 + p2);

        diagLen = hypot(p2(1)-p1(1), p2(2)-p1(2));
        radius_px = (1/(2*sqrt(2))) * diagLen;
        return;
    end
end

function [bestIdx, bestPts] = find_best_point_line(lines)
    bestIdx = nan; bestPts = []; bestCount = 0;
    for i = 1:numel(lines)
        s = char(lines(i));
        if ~contains(s,'(') || ~contains(s,',') || ~contains(s,')'), continue; end
        pts = parse_points_in_parentheses(s);
        if size(pts,1) > bestCount
            bestCount = size(pts,1);
            bestIdx = i;
            bestPts = pts;
        end
    end
    if bestCount < 2
        bestIdx = nan; bestPts = [];
    end
end

function pts = parse_points_in_parentheses(s)
    pts = [];
    tok = regexp(s, '\(\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*,\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*\)', 'tokens');
    if isempty(tok), return; end
    pts = zeros(numel(tok), 2);
    for i = 1:numel(tok)
        pts(i,1) = str2double(tok{i}{1});
        pts(i,2) = str2double(tok{i}{2});
    end
end

function [Xg, Yg, Ug, Vg] = vector_to_grid(x, y, u, v)
    ux = unique(x); uy = unique(y);
    ux = sort(ux);  uy = sort(uy);
    nx = numel(ux); ny = numel(uy);

    [~, ix] = ismember(x, ux);
    [~, iy] = ismember(y, uy);

    Ug = nan(ny, nx);
    Vg = nan(ny, nx);

    lin = sub2ind([ny,nx], iy, ix);
    Ug(lin) = u;
    Vg(lin) = v;

    [Xg, Yg] = meshgrid(ux, uy);
end
