%% ROI_Overlay_plus_NonAffine_D2min_CUMDISP_X_vs_nonX.m
% SAME as your affine (D2min) overlay/quant script, BUT:
%   - D2min is computed from CUMULATIVE displacement, not instantaneous.
%
% Meaning:
%   For each experiment, we read PIV frames in time order.
%   Each PIV file gives a displacement-per-frame field u(x,y), v(x,y) in um/frame.
%   We accumulate:
%       Ucum_k = sum_{j<=k} u_j
%       Vcum_k = sum_{j<=k} v_j
%   Then compute D2min on (Ucum_k, Vcum_k) at each k.
%
% Everything else (ROI parsing, overlay, output CSV) stays the same structure.

clearvars; close all; clc;

%% ================== USER SETTINGS ==================
RUN_MODE = 'nonX';  % 'nonX' | 'X' | 'both'

Experiments_nonX = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,22,24,25];
Experiments_X    = [3 7 8 9 10 14 15 16 17 24 25 26 27 28 29 30 31 32 33 34 35];

BASE_DIR_X    = 'PATH_TO_DATA/Experiments/SingleLightSpot_Data&PIV/';
BASE_DIR_nonX = 'PATH_TO_DATA/Experiments/NonX/';

% ---- PIV file convention ----
PIV_STEP = 2;              % PIV_ti - t(i+2)
POS_COL_X_PX = 1;          % x position (pixels) in PIV txt
POS_COL_Y_PX = 2;          % y position (pixels)
U_COL = 13;                % ux0 column
V_COL = 14;                % uy0 column

% If x/y units parsing fails:
DEFAULT_X_UNIT = 0.1;      % um/px
DEFAULT_Y_UNIT = 0.1;      % um/px

% ---- nonX Excel ----
XLSX_PATH  = 'ExperimentInfo.xlsx';
XLSX_SHEET = 1;
COL_ID     = 1;            % 'set3'
COL_CX_PX  = 3;            % CenterX(px)
COL_CY_PX  = 4;            % CenterY(px)
COL_R_UM   = 5;            % Radius(um)  (already um)

% ---- Plot settings ----
QUIVER_STRIDE = 2;         % plot every Nth vector in each direction
QUIVER_SCALE  = 0;         % 0 = autoscale
VEC_COLOR     = 'y';
ROI_COLOR     = 'c';

% ---- D2min settings ----
NBHW          = 3;         % neighborhood half-width: 1 => 3x3, 2 => 5x5
NORMALIZE_D2  = true;      % D2 / mean(|dr|^2) in neighborhood (dimensionless)
COMPUTE_ALL_FRAMES_D2 = true;   % if true: compute cumulative-D2 at ALL frames; else: only last

% ---- Visualize D2min on last PIV frame (heatmap) ----
SHOW_D2_HEATMAP_LAST = true;    % overlays D2 heatmap (semi-transparent)
D2_HEATMAP_ALPHA     = 0.35;    % transparency for heatmap overlay

SAVE_DIR_ROOT = fullfile(pwd, 'ROI_PIV_Overlays_plusD2_CumDisp');
DEBUG_PARSE = false;

%% ================== PREP ==================
if ~exist(SAVE_DIR_ROOT, 'dir'), mkdir(SAVE_DIR_ROOT); end

% Load Excel only if we will run nonX
T = table(); idCol = strings(0,1);
if any(strcmpi(RUN_MODE, {'nonx','both'}))
    if ~isfile(XLSX_PATH)
        error('RUN_MODE includes nonX but Excel not found: %s', XLSX_PATH);
    end
    T = readtable(XLSX_PATH, 'Sheet', XLSX_SHEET, 'VariableNamingRule','preserve');
    idCol = string(T{:, COL_ID});
end

% results table
Results = table();
Results.Group = strings(0,1);
Results.SetID = zeros(0,1);
Results.LastPIV_i = nan(0,1);
Results.D2_ROImean_last = nan(0,1);
Results.D2_ROImean_max = nan(0,1);
Results.D2_ROImean_mean = nan(0,1);
Results.ValidFrames = nan(0,1);

%% ================== RUN ==================
switch lower(RUN_MODE)
    case 'nonx'
        Results = [Results; run_one_group( ...
            'nonX', Experiments_nonX, BASE_DIR_nonX, true, ...
            SAVE_DIR_ROOT, ...
            T, idCol, ...
            COL_CX_PX, COL_CY_PX, COL_R_UM, ...
            PIV_STEP, POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
            DEFAULT_X_UNIT, DEFAULT_Y_UNIT, ...
            QUIVER_STRIDE, QUIVER_SCALE, VEC_COLOR, ROI_COLOR, ...
            NBHW, NORMALIZE_D2, COMPUTE_ALL_FRAMES_D2, ...
            SHOW_D2_HEATMAP_LAST, D2_HEATMAP_ALPHA, ...
            DEBUG_PARSE)];

    case 'x'
        Results = [Results; run_one_group( ...
            'X', Experiments_X, BASE_DIR_X, false, ...
            SAVE_DIR_ROOT, ...
            T, idCol, ...
            COL_CX_PX, COL_CY_PX, COL_R_UM, ...
            PIV_STEP, POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
            DEFAULT_X_UNIT, DEFAULT_Y_UNIT, ...
            QUIVER_STRIDE, QUIVER_SCALE, VEC_COLOR, ROI_COLOR, ...
            NBHW, NORMALIZE_D2, COMPUTE_ALL_FRAMES_D2, ...
            SHOW_D2_HEATMAP_LAST, D2_HEATMAP_ALPHA, ...
            DEBUG_PARSE)];

    case 'both'
        Results = [Results; run_one_group( ...
            'nonX', Experiments_nonX, BASE_DIR_nonX, true, ...
            SAVE_DIR_ROOT, ...
            T, idCol, ...
            COL_CX_PX, COL_CY_PX, COL_R_UM, ...
            PIV_STEP, POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
            DEFAULT_X_UNIT, DEFAULT_Y_UNIT, ...
            QUIVER_STRIDE, QUIVER_SCALE, VEC_COLOR, ROI_COLOR, ...
            NBHW, NORMALIZE_D2, COMPUTE_ALL_FRAMES_D2, ...
            SHOW_D2_HEATMAP_LAST, D2_HEATMAP_ALPHA, ...
            DEBUG_PARSE)];

        Results = [Results; run_one_group( ...
            'X', Experiments_X, BASE_DIR_X, false, ...
            SAVE_DIR_ROOT, ...
            T, idCol, ...
            COL_CX_PX, COL_CY_PX, COL_R_UM, ...
            PIV_STEP, POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
            DEFAULT_X_UNIT, DEFAULT_Y_UNIT, ...
            QUIVER_STRIDE, QUIVER_SCALE, VEC_COLOR, ROI_COLOR, ...
            NBHW, NORMALIZE_D2, COMPUTE_ALL_FRAMES_D2, ...
            SHOW_D2_HEATMAP_LAST, D2_HEATMAP_ALPHA, ...
            DEBUG_PARSE)];
    otherwise
        error('RUN_MODE must be ''nonX'', ''X'', or ''both''.');
end

% save results
outCSV = fullfile(SAVE_DIR_ROOT, sprintf('D2_results_%s_cumdisp.csv', RUN_MODE));
writetable(Results, outCSV);
fprintf('Saved results CSV: %s\n', outCSV);
fprintf('Done. Overlays saved under: %s\n', SAVE_DIR_ROOT);

%% ================== FUNCTIONS ==================

function Res = run_one_group(label, expList, baseDir, useExcel, ...
                             SAVE_DIR_ROOT, ...
                             T, idCol, ...
                             COL_CX_PX, COL_CY_PX, COL_R_UM, ...
                             PIV_STEP, POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
                             DEFAULT_X_UNIT, DEFAULT_Y_UNIT, ...
                             QUIVER_STRIDE, QUIVER_SCALE, VEC_COLOR, ROI_COLOR, ...
                             NBHW, NORMALIZE_D2, COMPUTE_ALL_FRAMES_D2, ...
                             SHOW_D2_HEATMAP_LAST, D2_HEATMAP_ALPHA, ...
                             DEBUG_PARSE)

    outDir = fullfile(SAVE_DIR_ROOT, label);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    Res = table();
    Res.Group = strings(0,1);
    Res.SetID = zeros(0,1);
    Res.LastPIV_i = nan(0,1);
    Res.D2_ROImean_last = nan(0,1);
    Res.D2_ROImean_max = nan(0,1);
    Res.D2_ROImean_mean = nan(0,1);
    Res.ValidFrames = nan(0,1);

    for ii = 1:numel(expList)
        setID = expList(ii);
        expDir = fullfile(baseDir, sprintf('set%d', setID));
        pivDir = fullfile(expDir, 'PIV_Results');

        if ~isfolder(expDir)
            warning('%s set%d: missing folder. Skipping.', label, setID);
            continue;
        end
        if ~isfolder(pivDir)
            warning('%s set%d: missing PIV_Results. Skipping.', label, setID);
            continue;
        end

        % parse x_unit/y_unit
        [x_unit, y_unit] = parse_xy_units_from_any_txt(expDir);
        if isnan(x_unit), x_unit = DEFAULT_X_UNIT; end
        if isnan(y_unit), y_unit = DEFAULT_Y_UNIT; end

        % load last TIFF
        tiffPath = find_first_tiff(expDir);
        if isempty(tiffPath)
            warning('%s set%d: no TIFF found. Skipping.', label, setID);
            continue;
        end
        info = imfinfo(tiffPath);
        lastFrame = numel(info);
        I = imread(tiffPath, 'Index', lastFrame);

        % ROI in pixels (for overlay) + ROI in um (for D2 masking on um-grid)
        if useExcel
            key = "set" + string(setID);
            row = find(idCol == key, 1, 'first');
            if isempty(row)
                warning('%s set%d: not in Excel. Skipping.', label, setID);
                continue;
            end
            cx_px = T{row, COL_CX_PX};
            cy_px = T{row, COL_CY_PX};
            r_um  = T{row, COL_R_UM};

            if any(isnan([cx_px, cy_px, r_um]))
                warning('%s set%d: Excel has NaNs. Skipping.', label, setID);
                continue;
            end

            center_px = [cx_px, cy_px];
            rx_px = r_um / x_unit;    % ellipse radii in pixels
            ry_px = r_um / y_unit;

            center_um = [cx_px * x_unit, cy_px * y_unit];
            radius_um = r_um;
        else
            [center_px, radius_px] = parse_roi_from_txt_pixels(expDir, DEBUG_PARSE, setID);
            if any(isnan([center_px, radius_px]))
                warning('%s set%d: failed ROI parse. Skipping.', label, setID);
                continue;
            end
            rx_px = radius_px; ry_px = radius_px;

            center_um = [center_px(1) * x_unit, center_px(2) * y_unit];
            radius_um = radius_px * sqrt(x_unit * y_unit); % approx
        end

        % list PIV i in order
        i_list = list_all_piv_i(pivDir, PIV_STEP);
        if isempty(i_list)
            warning('%s set%d: no PIV files. Will save ROI-only overlay.', label, setID);
            last_i = nan;
        else
            last_i = i_list(end);
        end

        % --- compute D2 from CUMULATIVE displacement ---
        D2_last_roiMean = nan;
        D2_roiMean_max  = nan;
        D2_roiMean_mean = nan;
        validFrames     = 0;

        % for heatmap overlay
        D2_last = []; Xum_last = []; Yum_last = [];

        if ~isnan(last_i)
            if COMPUTE_ALL_FRAMES_D2
                [roiMeans, D2_last, Xum_last, Yum_last] = compute_D2_cumdisp_over_time( ...
                    pivDir, i_list, PIV_STEP, ...
                    POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
                    x_unit, y_unit, ...
                    center_um, radius_um, ...
                    NBHW, NORMALIZE_D2);

                roiMeans = roiMeans(isfinite(roiMeans));
                validFrames = numel(roiMeans);
                if validFrames > 0
                    D2_last_roiMean = roiMeans(end);
                    D2_roiMean_max  = max(roiMeans);
                    D2_roiMean_mean = mean(roiMeans);
                end
            else
                % only compute LAST cumulative state (still accumulates internally)
                [D2_last, Xum_last, Yum_last, ok_last, D2_last_roiMean] = compute_D2_cumdisp_last( ...
                    pivDir, i_list, PIV_STEP, ...
                    POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
                    x_unit, y_unit, ...
                    center_um, radius_um, ...
                    NBHW, NORMALIZE_D2);
                if ok_last
                    D2_roiMean_max  = D2_last_roiMean;
                    D2_roiMean_mean = D2_last_roiMean;
                    validFrames     = 1;
                end
            end
        end

        % --- make overlay figure ---
        fig = figure('Color','w', 'Visible','off');
        imshow(I, [], 'InitialMagnification', 'fit'); hold on;

        if isnan(last_i)
            title(sprintf('%s set%d | TIFF %d | PIV: none', label, setID, lastFrame));
        else
            title(sprintf('%s set%d | TIFF %d | last PIV i=%d | D2cum_last=%.3g | D2cum_max=%.3g', ...
                label, setID, lastFrame, last_i, D2_last_roiMean, D2_roiMean_max));
        end

        % ROI overlay
        if abs(rx_px - ry_px) > 1e-3*max(rx_px, ry_px)
            plot_ellipse(center_px(1), center_px(2), rx_px, ry_px, ROI_COLOR, 2);
        else
            draw_circle(center_px(1), center_px(2), mean([rx_px,ry_px]), ROI_COLOR, 2);
        end

        % Overlay vectors (LAST instantaneous PIV) in pixel coords (visual check)
        if ~isnan(last_i)
            pivFile = fullfile(pivDir, sprintf('PIV_t%d-t%d.txt', last_i, last_i+PIV_STEP));
            if exist(pivFile, 'file') == 2
                data = load(pivFile);
                needCols = max([POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL]);
                if size(data,2) >= needCols
                    x_px = data(:, POS_COL_X_PX);
                    y_px = data(:, POS_COL_Y_PX);
                    u_px = data(:, U_COL) / PIV_STEP;
                    v_px = data(:, V_COL) / PIV_STEP;
                    [Xg, Yg, Ug, Vg] = vector_to_grid(x_px, y_px, u_px, v_px);

                    Xs = Xg(1:QUIVER_STRIDE:end, 1:QUIVER_STRIDE:end);
                    Ys = Yg(1:QUIVER_STRIDE:end, 1:QUIVER_STRIDE:end);
                    Us = Ug(1:QUIVER_STRIDE:end, 1:QUIVER_STRIDE:end);
                    Vs = Vg(1:QUIVER_STRIDE:end, 1:QUIVER_STRIDE:end);

                    quiver(Xs, Ys, Us, Vs, QUIVER_SCALE, VEC_COLOR, 'LineWidth', 1);
                end
            end
        end

        % Optional: D2 heatmap overlay for LAST cumulative state
        if SHOW_D2_HEATMAP_LAST && ~isnan(last_i) && ~isempty(D2_last)
            Xpx = Xum_last / x_unit;
            Ypx = Yum_last / y_unit;
            imagesc(Xpx(1,:), Ypx(:,1), D2_last);
            set(gca, 'YDir','reverse');
            colormap(gca, parula);
            alpha(D2_HEATMAP_ALPHA);
            uistack(findobj(gca,'Type','line'), 'top');
            uistack(findobj(gca,'Type','quiver'), 'top');
        end

        outPath = fullfile(outDir, sprintf('set%d_%s_overlay_cumdisp.png', setID, label));
        exportgraphics(fig, outPath, 'Resolution', 200);
        close(fig);
        fprintf('Saved: %s\n', outPath);

        Res = [Res; {string(label), setID, last_i, D2_last_roiMean, D2_roiMean_max, D2_roiMean_mean, validFrames}]; %#ok<AGROW>
        Res.Properties.VariableNames = {'Group','SetID','LastPIV_i','D2_ROImean_last','D2_ROImean_max','D2_ROImean_mean','ValidFrames'};
    end
end

function [roiMeans, D2_last, Xum_last, Yum_last] = compute_D2_cumdisp_over_time( ...
            pivDir, i_list, pivStep, POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
            x_unit, y_unit, center_um, radius_um, nbhw, normalizeD2)

    roiMeans = nan(size(i_list));

    Ucum = []; Vcum = [];
    roiMask = [];
    gridSet = false;

    D2_last = []; Xum_last = []; Yum_last = [];

    for kk = 1:numel(i_list)
        i = i_list(kk);

        [Xum, Yum, Ustep, Vstep, ok] = load_uv_grid_um_perframe( ...
            pivDir, i, pivStep, POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, x_unit, y_unit);
        if ~ok, continue; end

        if ~gridSet
            Ucum = zeros(size(Ustep), 'like', Ustep);
            Vcum = zeros(size(Vstep), 'like', Vstep);
            roiMask = hypot(Xum - center_um(1), Yum - center_um(2)) <= radius_um;
            gridSet = true;
        else
            if ~isequal(size(Ustep), size(Ucum))
                continue;
            end
        end

        Ucum = Ucum + Ustep;
        Vcum = Vcum + Vstep;

        D2 = nonaffinity_D2min_grid(Xum, Yum, Ucum, Vcum, nbhw, normalizeD2);

        vals = D2(roiMask);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            roiMeans(kk) = mean(vals);
        end

        if kk == numel(i_list)
            D2_last = D2;
            Xum_last = Xum;
            Yum_last = Yum;
        end
    end
end

function [D2_last, Xum_last, Yum_last, ok_last, roiMean_last] = compute_D2_cumdisp_last( ...
            pivDir, i_list, pivStep, POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
            x_unit, y_unit, center_um, radius_um, nbhw, normalizeD2)

    ok_last = false;
    roiMean_last = nan;

    Ucum = []; Vcum = [];
    roiMask = [];
    gridSet = false;

    D2_last = []; Xum_last = []; Yum_last = [];

    for kk = 1:numel(i_list)
        i = i_list(kk);

        [Xum, Yum, Ustep, Vstep, ok] = load_uv_grid_um_perframe( ...
            pivDir, i, pivStep, POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, x_unit, y_unit);
        if ~ok, continue; end

        if ~gridSet
            Ucum = zeros(size(Ustep), 'like', Ustep);
            Vcum = zeros(size(Vstep), 'like', Vstep);
            roiMask = hypot(Xum - center_um(1), Yum - center_um(2)) <= radius_um;
            gridSet = true;
        else
            if ~isequal(size(Ustep), size(Ucum))
                continue;
            end
        end

        Ucum = Ucum + Ustep;
        Vcum = Vcum + Vstep;

        if kk == numel(i_list)
            D2_last = nonaffinity_D2min_grid(Xum, Yum, Ucum, Vcum, nbhw, normalizeD2);
            Xum_last = Xum; Yum_last = Yum;

            vals = D2_last(roiMask);
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                roiMean_last = mean(vals);
                ok_last = true;
            end
        end
    end
end

function [Xum, Yum, Ug, Vg, ok] = load_uv_grid_um_perframe(pivDir, i, pivStep, ...
                                                          POS_COL_X_PX, POS_COL_Y_PX, U_COL, V_COL, ...
                                                          x_unit, y_unit)
    % Load ONE PIV file as displacement-per-frame in um/frame.
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

function tiffPath = find_first_tiff(expDir)
    tiffPath = '';
    tifFiles = dir(fullfile(expDir, '*.tif'));
    tifFiles = tifFiles(~startsWith({tifFiles.name}, '._'));
    if isempty(tifFiles), return; end
    tiffPath = fullfile(expDir, tifFiles(1).name);
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

function D2 = nonaffinity_D2min_grid(Xg, Yg, Ug, Vg, nbhw, normalizeD2)
    % Falk–Langer D2min on a grid. Here Ug,Vg can be instantaneous or cumulative.
    [ny, nx] = size(Ug);
    D2 = nan(ny, nx);

    for iy = 1:ny
        for ix = 1:nx
            ui = Ug(iy,ix); vi = Vg(iy,ix);
            if ~isfinite(ui) || ~isfinite(vi), continue; end

            y1 = max(1, iy-nbhw); y2 = min(ny, iy+nbhw);
            x1 = max(1, ix-nbhw); x2 = min(nx, ix+nbhw);

            xi = Xg(iy,ix);
            yi0 = Yg(iy,ix);

            DR = []; DV = [];

            for jy = y1:y2
                for jx = x1:x2
                    if jy==iy && jx==ix, continue; end
                    uj = Ug(jy,jx); vj = Vg(jy,jx);
                    if ~isfinite(uj) || ~isfinite(vj), continue; end

                    dx = Xg(jy,jx) - xi;
                    dy = Yg(jy,jx) - yi0;

                    du = uj - ui;
                    dv = vj - vi;

                    DR(:,end+1) = [dx; dy]; %#ok<AGROW>
                    DV(:,end+1) = [du; dv]; %#ok<AGROW>
                end
            end

            if size(DR,2) < 3, continue; end

            A = DR*DR.';
            if rcond(A) < 1e-12, continue; end

            E = (DV*DR.') / A;
            R = DV - E*DR;

            res2 = sum(R.^2, 1);
            d2 = mean(res2);

            if normalizeD2
                dr2 = sum(DR.^2, 1);
                d2 = d2 / max(mean(dr2), eps);
            end

            D2(iy,ix) = d2;
        end
    end
end

function draw_circle(cx, cy, r, colorChar, lw)
    th = linspace(0, 2*pi, 300);
    plot(cx + r*cos(th), cy + r*sin(th), colorChar, 'LineWidth', lw);
    plot(cx, cy, '+', 'Color', colorChar, 'LineWidth', lw);
end

function plot_ellipse(cx, cy, rx, ry, colorChar, lw)
    th = linspace(0, 2*pi, 300);
    plot(cx + rx*cos(th), cy + ry*sin(th), colorChar, 'LineWidth', lw);
    plot(cx, cy, '+', 'Color', colorChar, 'LineWidth', lw);
end
