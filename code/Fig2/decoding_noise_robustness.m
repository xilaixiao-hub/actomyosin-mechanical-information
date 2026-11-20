clearvars;
close all;
clc;

%% ============================================================
% Multi-spot pipeline-level noise injection analysis
%
% This script:
% 1. Imports multi-spot experiment image stack
% 2. Computes raw temporal STD map from the movie
% 3. Estimates raw noise scale from inactive bit ROIs
% 4. Adds synthetic raw noise to the raw STD map BEFORE normalization
% 5. Re-normalizes the noisy STD map
% 6. Extracts 32 site-level signals S_i
% 7. Decodes the noisy pattern
% 8. Calculates exact decoding ratio, bitwise accuracy, and recovered active sites
%
% Main idea:
%   UnINTENSITY_Std_noisy = UnINTENSITY_Std + alpha * noiseMap
%
% where noiseMap is newly generated for every trial.
%% ============================================================

%% ---------------- USER SETTINGS ----------------

Experiments = 1;

directory = ['../../Experiments/MultiLightSpots_Data&PIV/set' num2str(Experiments) '/'];

% True activated sites in the real experiment
answer = [5 6 7 8 13 17 19 23 24 29 30 32];

% Candidate positions used in your original decoding code
Inward_Index = [5 6 7 8 13 17 19 23 24 29 30 31 32];

% Decoder mode:
% 'KL'   = use original KL-style decoding
% 'TopK' = choose top numActive sites by noisy S_i
decoderMode = 'KL';

% Noise settings
numTrials = 2000;          % increase to 5000 later if runtime is okay
noiseFolds = 0:0.25:5;     % alpha = additional raw noise fold

% Noise type:
% 'pixel'  = independent pixel-level noise
% 'smooth' = spatially correlated noise map
noiseType = 'smooth';

% Smoothing length for smooth noise, in pixels
smoothSigma_px = 3;

% Random seed
rng(1);

% Output folder
outDir = ['PipelineNoiseInjection_set' num2str(Experiments) '_' decoderMode '_' noiseType];
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ---------------- READ METADATA TXT FILE ----------------

txtFiles = dir(fullfile(directory, '*.txt'));
txtFiles = txtFiles(~startsWith({txtFiles.name}, '._'));

if isempty(txtFiles)
    error('No txt metadata file found in directory: %s', directory);
end

filePath_txt = fullfile(directory, txtFiles(1).name);

fileID = fopen(filePath_txt, 'r');
C = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);

lines = C{1,1};

% Find x and y units
line2 = lines{5, 1};
parts = strsplit(line2, ':');
start_str = strtrim(parts{2});
step_str = strsplit(start_str, '*');
x_unit = str2double(strtrim(step_str{2}));

line3 = lines{6, 1};
parts = strsplit(line3, ':');
start_str = strtrim(parts{2});
step_str = strsplit(start_str, '*');
y_unit = str2double(strtrim(step_str{2}));

fprintf('x unit = %.5f um/pixel\n', x_unit);
fprintf('y unit = %.5f um/pixel\n', y_unit);

%% ---------------- FIND LIGHT RADIUS FROM METADATA ----------------

CENTERS = zeros(12,2);
LRadius = [];

for p = 1:12

    line = lines{end-3-p, 1};

    start_idx = strfind(line, '(');
    end_idx = strfind(line, ')');

    points = [];

    for ii = 1:length(start_idx)

        point_str = line(start_idx(ii)+1:end_idx(ii)-1);
        coords = strsplit(point_str, ',');

        if length(coords) == 2
            point = [str2double(strtrim(coords{1})), ...
                     str2double(strtrim(coords{2}))];
            points = [points; point];
        end
    end

    center = [0.5*(points(1,1)+points(2,1))*x_unit, ...
              0.5*(points(1,2)+points(2,2))*y_unit];

    CENTERS(p,:) = center;

    LRadius(end+1) = (1/(2*sqrt(2))) * sqrt( ...
        ((points(2,1)-points(1,1))*x_unit)^2 + ...
        ((points(2,2)-points(1,2))*y_unit)^2 );
end

R_light = mean(LRadius);
R_light = round(R_light * 1000) / 1000;

fprintf('Mean light radius = %.3f um\n', R_light);

%% ---------------- IMPORT TIFF STACK ----------------

tifFiles = dir(fullfile(directory, '*.tif'));
tifFiles = tifFiles(~startsWith({tifFiles.name}, '._'));

if isempty(tifFiles)
    error('No tif file found in directory: %s', directory);
end

filePath_tif = fullfile(directory, tifFiles(1).name);

info = imfinfo(filePath_tif);
numFrames = numel(info);

fprintf('Using TIFF file: %s\n', tifFiles(1).name);
fprintf('Number of TIFF frames = %d\n', numFrames);

I = imread(filePath_tif, 'Index', 1);

unique_x = (1:size(I,2)) * x_unit;
unique_y = (1:size(I,1)) * y_unit;

[X, Y] = meshgrid(unique_x, unique_y);

Intensity_Array = zeros(size(I,1), size(I,2), numFrames, 'uint16');

for i = 1:numFrames
    Intensity_Array(:,:,i) = imread(filePath_tif, 'Index', i);
end

%% ---------------- COMPUTE RAW STD MAP ----------------

% Raw temporal STD map, before any normalization
UnINTENSITY_Std = std(double(Intensity_Array), 0, 3);

% Original normalized STD map, same as your original pipeline
OriginalNormSTD = normalize_with_threshold(UnINTENSITY_Std);

fprintf('Raw STD map calculated.\n');

%% ---------------- DEFINE BIT GRID ----------------

Row_coords = [53.9 73 92.45 111.9 133 152.25 169.875 187.5];
Colum_coords = [65 87 108 128];

N = 32;

B_true = zeros(N,1);
B_true(answer) = 1;

numActive = sum(B_true);

fprintf('True active sites: ');
fprintf('%d ', answer);
fprintf('\n');
fprintf('Number of active sites = %d\n', numActive);

%% ---------------- EXTRACT ORIGINAL SITE SIGNALS ----------------

% Raw site values from raw STD map
S_raw_original = extract_site_values_from_map( ...
    UnINTENSITY_Std, X, Y, Row_coords, Colum_coords, R_light);

S_raw_original = S_raw_original(:);

% Normalized site values from original normalized STD map
S_original = extract_site_values_from_map( ...
    OriginalNormSTD, X, Y, Row_coords, Colum_coords, R_light);

S_original = S_original(:);

% Raw noise scale estimated from inactive sites in raw STD units
sigma_raw_noise = std(S_raw_original(B_true == 0), 'omitnan');

fprintf('Raw inactive-site sigma_noise = %.5f\n', sigma_raw_noise);

Delta_original = mean(S_original(B_true == 1), 'omitnan') - ...
                 mean(S_original(B_true == 0), 'omitnan');

fprintf('Original normalized Delta = %.5f\n', Delta_original);

%% ---------------- BUILD CANDIDATE CONFIGURATIONS FOR KL DECODER ----------------

if strcmpi(decoderMode, 'KL')

    candidateConfigs = {};
    candidateMatrices = {};

    nCandPos = length(Inward_Index);

    for k = 1:nCandPos

        index_combinations = nchoosek(1:nCandPos, k);

        for oo = 1:size(index_combinations,1)

            candidateIdx = Inward_Index(index_combinations(oo,:));

            Possible_STD = zeros(1,N);
            Possible_STD(candidateIdx) = 1;

            Possible_STD = reshape(Possible_STD, 4, 8);
            Possible_STD = Possible_STD';

            Possible_STD = normalize_probability_matrix(Possible_STD);

            candidateConfigs{end+1,1} = candidateIdx;
            candidateMatrices{end+1,1} = Possible_STD;
        end
    end

    fprintf('Number of KL candidate configurations = %d\n', length(candidateConfigs));
end

%% ---------------- CHECK ORIGINAL UNPERTURBED DECODING ----------------

if strcmpi(decoderMode, 'KL')

    decoded_original = decode_KL_from_signal(S_original, candidateConfigs, candidateMatrices);

elseif strcmpi(decoderMode, 'TopK')

    decoded_original = decode_TopK_from_signal(S_original, numActive);

else
    error('Unknown decoderMode. Use KL or TopK.');
end

B_decoded_original = zeros(N,1);
B_decoded_original(decoded_original) = 1;

exact_original = all(B_decoded_original == B_true);
bitwise_original = mean(B_decoded_original == B_true);
recovered_original = sum(B_decoded_original == 1 & B_true == 1);

fprintf('\nOriginal unperturbed decoding:\n');
fprintf('Decoded sites: ');
fprintf('%d ', decoded_original);
fprintf('\n');
fprintf('Exact correct = %d\n', exact_original);
fprintf('Bitwise accuracy = %.3f\n', bitwise_original);
fprintf('Recovered active sites = %d/%d\n', recovered_original, numActive);

%% ---------------- SAVE ORIGINAL STD FIGURES ----------------

figure('Color','w');
imagesc(unique_x, unique_y, OriginalNormSTD);
axis image;
set(gca, 'YDir', 'normal');
colorbar;
title('Original normalized STD map');
xlabel('X (\mum)');
ylabel('Y (\mum)');
hold on;
plot_bit_circles(Row_coords, Colum_coords, R_light, answer, []);
saveas(gcf, fullfile(outDir, 'Original_Normalized_STD_Map.png'));

figure('Color','w');
imagesc(unique_x, unique_y, UnINTENSITY_Std);
axis image;
set(gca, 'YDir', 'normal');
colorbar;
title('Original raw STD map');
xlabel('X (\mum)');
ylabel('Y (\mum)');
hold on;
plot_bit_circles(Row_coords, Colum_coords, R_light, answer, []);
saveas(gcf, fullfile(outDir, 'Original_Raw_STD_Map.png'));

%% ---------------- PIPELINE-LEVEL NOISE INJECTION LOOP ----------------

fprintf('\n===== Running raw-map noise injection before normalization =====\n');
fprintf('Noise type: %s\n', noiseType);

exactCorrectRatio = nan(size(noiseFolds));
bitwiseCorrectRatio = nan(size(noiseFolds));
meanNumCorrectActive = nan(size(noiseFolds));
stdNumCorrectActive = nan(size(noiseFolds));
meanDeltaNoisy = nan(size(noiseFolds));
stdDeltaNoisy = nan(size(noiseFolds));

% Store example noisy maps for alpha = 1 if present
saveExampleAlpha = 1;
exampleSaved = false;

for aa = 1:numel(noiseFolds)

    alpha = noiseFolds(aa);

    exactCorrect = false(numTrials,1);
    bitwiseAccuracy = nan(numTrials,1);
    numCorrectActive = nan(numTrials,1);
    DeltaNoisy = nan(numTrials,1);

    for tt = 1:numTrials

        % Generate new raw noise map for every trial
        noiseMap = generate_noise_map(size(UnINTENSITY_Std), sigma_raw_noise, noiseType, smoothSigma_px);

        % Add alpha-fold noise to raw STD map BEFORE normalization
        NoisyRawSTD = UnINTENSITY_Std + alpha * noiseMap;

        % STD cannot be negative
        NoisyRawSTD(NoisyRawSTD < 0) = 0;

        % Normalize noisy raw STD map using the same original function
        NoisyNormSTD = normalize_with_threshold(NoisyRawSTD);

        % Extract 32 site-level signals from noisy normalized STD map
        S_noisy = extract_site_values_from_map( ...
            NoisyNormSTD, X, Y, Row_coords, Colum_coords, R_light);

        S_noisy = S_noisy(:);

        % Decode noisy site-level signal
        if strcmpi(decoderMode, 'KL')

            decodedIdx = decode_KL_from_signal(S_noisy, candidateConfigs, candidateMatrices);

        elseif strcmpi(decoderMode, 'TopK')

            decodedIdx = decode_TopK_from_signal(S_noisy, numActive);

        end

        B_decoded = zeros(N,1);
        B_decoded(decodedIdx) = 1;

        exactCorrect(tt) = all(B_decoded == B_true);
        bitwiseAccuracy(tt) = mean(B_decoded == B_true);
        numCorrectActive(tt) = sum(B_decoded == 1 & B_true == 1);

        DeltaNoisy(tt) = mean(S_noisy(B_true == 1), 'omitnan') - ...
                         mean(S_noisy(B_true == 0), 'omitnan');

        % Save one example noisy map near alpha = 1
        if ~exampleSaved && abs(alpha - saveExampleAlpha) < 1e-12 && tt == 1

            figure('Color','w');
            imagesc(unique_x, unique_y, NoisyRawSTD);
            axis image;
            set(gca, 'YDir', 'normal');
            colorbar;
            title(['Example noisy raw STD map, \alpha = ' num2str(alpha)]);
            xlabel('X (\mum)');
            ylabel('Y (\mum)');
            hold on;
            plot_bit_circles(Row_coords, Colum_coords, R_light, answer, []);
            saveas(gcf, fullfile(outDir, ['Example_NoisyRawSTD_alpha' num2str(alpha) '.png']));

            figure('Color','w');
            imagesc(unique_x, unique_y, NoisyNormSTD);
            axis image;
            set(gca, 'YDir', 'normal');
            colorbar;
            title(['Example noisy normalized STD map, \alpha = ' num2str(alpha)]);
            xlabel('X (\mum)');
            ylabel('Y (\mum)');
            hold on;
            plot_bit_circles(Row_coords, Colum_coords, R_light, answer, []);
            saveas(gcf, fullfile(outDir, ['Example_NoisyNormSTD_alpha' num2str(alpha) '.png']));

            exampleSaved = true;
        end
    end

    exactCorrectRatio(aa) = mean(exactCorrect);
    bitwiseCorrectRatio(aa) = mean(bitwiseAccuracy, 'omitnan');
    meanNumCorrectActive(aa) = mean(numCorrectActive, 'omitnan');
    stdNumCorrectActive(aa) = std(numCorrectActive, 'omitnan');
    meanDeltaNoisy(aa) = mean(DeltaNoisy, 'omitnan');
    stdDeltaNoisy(aa) = std(DeltaNoisy, 'omitnan');

    fprintf('Noise fold %.2f: exact = %.3f, bitwise = %.3f, recovered active = %.2f/%d, Delta = %.4f\n', ...
        alpha, exactCorrectRatio(aa), bitwiseCorrectRatio(aa), ...
        meanNumCorrectActive(aa), numActive, meanDeltaNoisy(aa));
end

%% ---------------- TOTAL EFFECTIVE NOISE FOLD ----------------

% alpha is added independent noise relative to raw inactive-site variability.
% Since the original data already contains baseline noise, a rough effective
% total noise fold is sqrt(1 + alpha^2).
totalNoiseFold = sqrt(1 + noiseFolds.^2);

%% ---------------- SAVE RESULT TABLE ----------------

ResultTable = table( ...
    noiseFolds(:), ...
    totalNoiseFold(:), ...
    exactCorrectRatio(:), ...
    bitwiseCorrectRatio(:), ...
    meanNumCorrectActive(:), ...
    stdNumCorrectActive(:), ...
    meanDeltaNoisy(:), ...
    stdDeltaNoisy(:), ...
    'VariableNames', { ...
        'AdditionalNoiseFold_alpha', ...
        'EstimatedTotalNoiseFold', ...
        'ExactCorrectRatio', ...
        'BitwiseCorrectRatio', ...
        'MeanRecoveredActiveSites', ...
        'StdRecoveredActiveSites', ...
        'MeanDeltaNoisy', ...
        'StdDeltaNoisy' ...
    });

outTableName = fullfile(outDir, ['PipelineNoiseInjection_Result_set' num2str(Experiments) '_' decoderMode '_' noiseType '.csv']);
writetable(ResultTable, outTableName);

fprintf('\nSaved result table:\n%s\n', outTableName);

%% ---------------- PLOTS: X AXIS = ADDITIONAL NOISE FOLD ----------------

FZ = 16;

figure('Color','w');
plot(noiseFolds, exactCorrectRatio, '-o', ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7);
xlabel('Additional raw noise fold before normalization, \alpha', 'FontSize', FZ);
ylabel('Exact correct decoding ratio', 'FontSize', FZ);
title(['Pipeline-level noise injection, set' num2str(Experiments) ', ' decoderMode], 'FontSize', FZ);
ylim([-0.05 1.05]);
xlim([min(noiseFolds) max(noiseFolds)]);
set(gca, 'FontSize', FZ, 'LineWidth', 1.2);
box on;
saveas(gcf, fullfile(outDir, 'ExactCorrectRatio_vs_AdditionalNoiseFold.png'));

figure('Color','w');
plot(noiseFolds, bitwiseCorrectRatio, '-o', ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7);
xlabel('Additional raw noise fold before normalization, \alpha', 'FontSize', FZ);
ylabel('Mean bitwise correct ratio', 'FontSize', FZ);
title(['Bitwise robustness, set' num2str(Experiments) ', ' decoderMode], 'FontSize', FZ);
ylim([0.4 1.05]);
xlim([min(noiseFolds) max(noiseFolds)]);
set(gca, 'FontSize', FZ, 'LineWidth', 1.2);
box on;
saveas(gcf, fullfile(outDir, 'BitwiseCorrectRatio_vs_AdditionalNoiseFold.png'));

figure('Color','w');
errorbar(noiseFolds, meanNumCorrectActive, stdNumCorrectActive, '-o', ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7, ...
    'CapSize', 10);
xlabel('Additional raw noise fold before normalization, \alpha', 'FontSize', FZ);
ylabel('Recovered true active sites', 'FontSize', FZ);
title(['Recovered active sites, set' num2str(Experiments) ', ' decoderMode], 'FontSize', FZ);
ylim([0 numActive + 1]);
xlim([min(noiseFolds) max(noiseFolds)]);
set(gca, 'FontSize', FZ, 'LineWidth', 1.2);
box on;
saveas(gcf, fullfile(outDir, 'RecoveredActiveSites_vs_AdditionalNoiseFold.png'));

%% ---------------- PLOTS: X AXIS = ESTIMATED TOTAL NOISE FOLD ----------------

figure('Color','w');
plot(totalNoiseFold, exactCorrectRatio, '-o', ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7);
xlabel('Estimated total noise fold relative to original data', 'FontSize', FZ);
ylabel('Exact correct decoding ratio', 'FontSize', FZ);
title(['Pipeline-level robustness, set' num2str(Experiments) ', ' decoderMode], 'FontSize', FZ);
ylim([-0.05 1.05]);
xlim([min(totalNoiseFold) max(totalNoiseFold)]);
set(gca, 'FontSize', FZ, 'LineWidth', 1.2);
box on;
saveas(gcf, fullfile(outDir, 'ExactCorrectRatio_vs_TotalNoiseFold.png'));

figure('Color','w');
plot(totalNoiseFold, bitwiseCorrectRatio, '-o', ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7);
xlabel('Estimated total noise fold relative to original data', 'FontSize', FZ);
ylabel('Mean bitwise correct ratio', 'FontSize', FZ);
title(['Bitwise robustness, set' num2str(Experiments) ', ' decoderMode], 'FontSize', FZ);
ylim([0.4 1.05]);
xlim([min(totalNoiseFold) max(totalNoiseFold)]);
set(gca, 'FontSize', FZ, 'LineWidth', 1.2);
box on;
saveas(gcf, fullfile(outDir, 'BitwiseCorrectRatio_vs_TotalNoiseFold.png'));

figure('Color','w');
errorbar(totalNoiseFold, meanNumCorrectActive, stdNumCorrectActive, '-o', ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7, ...
    'CapSize', 10);
xlabel('Estimated total noise fold relative to original data', 'FontSize', FZ);
ylabel('Recovered true active sites', 'FontSize', FZ);
title(['Recovered active sites, set' num2str(Experiments) ', ' decoderMode], 'FontSize', FZ);
ylim([0 numActive + 1]);
xlim([min(totalNoiseFold) max(totalNoiseFold)]);
set(gca, 'FontSize', FZ, 'LineWidth', 1.2);
box on;
saveas(gcf, fullfile(outDir, 'RecoveredActiveSites_vs_TotalNoiseFold.png'));

fprintf('\nDone. Results saved in folder:\n%s\n', outDir);

%% ============================================================
% LOCAL HELPER FUNCTIONS
% Keep all helper functions at the very end of the script.
%% ============================================================

function normalized_matrix = normalize_with_threshold(matrix)

    matrix = double(matrix);

    pd = fitdist(matrix(:), 'Normal');

    mean_val = pd.mu;
    std_dev = pd.sigma;

    upper_threshold = mean_val + 3 * std_dev;
    lower_threshold = mean_val - 3 * std_dev;

    matrix(matrix > upper_threshold) = upper_threshold;
    matrix(matrix < lower_threshold) = lower_threshold;

    denom = max(matrix(:)) - min(matrix(:));

    if denom == 0 || isnan(denom)
        normalized_matrix = zeros(size(matrix));
    else
        normalized_matrix = (matrix - min(matrix(:))) ./ denom;
    end
end

function S = extract_site_values_from_map(STDmap, X, Y, Row_coords, Colum_coords, R_light)

    S = [];

    for m = 1:length(Row_coords)

        for n = 1:length(Colum_coords)

            center_x = Colum_coords(n);
            center_y = Row_coords(m);

            distance_from_center = sqrt((X - center_x).^2 + ...
                                        (Y - center_y).^2);

            bit_mask = distance_from_center <= R_light;

            values = STDmap(bit_mask);

            S(end+1) = mean(values, 'omitnan');
        end
    end
end

function noiseMap = generate_noise_map(mapSize, sigma_raw_noise, noiseType, smoothSigma_px)

    rawNoise = randn(mapSize);

    if strcmpi(noiseType, 'pixel')

        noiseMap = rawNoise;

    elseif strcmpi(noiseType, 'smooth')

        % Spatially correlated noise.
        % Requires Image Processing Toolbox for imgaussfilt.
        noiseMap = imgaussfilt(rawNoise, smoothSigma_px);

    else
        error('Unknown noiseType. Use pixel or smooth.');
    end

    % Normalize the generated noise map to have the desired STD
    currentStd = std(noiseMap(:), 'omitnan');

    if currentStd == 0 || isnan(currentStd)
        noiseMap = zeros(mapSize);
    else
        noiseMap = noiseMap / currentStd * sigma_raw_noise;
    end
end

function probMatrix = normalize_probability_matrix(M)

    M = double(M);

    M(M < 0) = 0;

    totalSum = sum(M(:));

    if totalSum <= 0 || isnan(totalSum)
        probMatrix = ones(size(M)) / numel(M);
    else
        probMatrix = M / totalSum;
    end
end

function D_KL = calculate_kl_divergence_local(P, Q)

    assert(all(size(P) == size(Q)), ...
        'Matrices P and Q must have the same size');

    P = double(P);
    Q = double(Q);

    P(P <= 0) = eps;
    Q(Q <= 0) = eps;

    P = P / sum(P(:));
    Q = Q / sum(Q(:));

    D_KL = sum(P .* log(P ./ Q), 'all');
end

function decodedIdx = decode_TopK_from_signal(S_signal, numActive)

    [~, sortedIdx] = sort(S_signal, 'descend');
    decodedIdx = sort(sortedIdx(1:numActive));
end

function decodedIdx = decode_KL_from_signal(S_signal, candidateConfigs, candidateMatrices)

    % Convert 32-site signal vector into same 8x4 matrix format
    % as your original decoding code.
    M = reshape(S_signal(:), 4, 8);
    M = M';

    % Match your original decoding preprocessing
    M = M.^2;

    M_prob = normalize_probability_matrix(M);

    KL_scores = nan(length(candidateMatrices),1);

    for cc = 1:length(candidateMatrices)

        Possible_STD = candidateMatrices{cc};

        % Same direction as your original code:
        % calculate_kl_divergence(Possible_STD, INTENSITY_STD)
        KL_scores(cc) = calculate_kl_divergence_local(Possible_STD, M_prob);
    end

    bestIdx = find(KL_scores == min(KL_scores), 1, 'first');

    decodedIdx = sort(candidateConfigs{bestIdx});
end

function plot_bit_circles(Row_coords, Colum_coords, R_light, activeIdx, inactiveIdx)

    if nargin < 5
        inactiveIdx = [];
    end

    theta = linspace(0, 2*pi, 100);

    running = 0;

    for m = 1:length(Row_coords)

        for n = 1:length(Colum_coords)

            running = running + 1;

            cx = Colum_coords(n);
            cy = Row_coords(m);

            xCircle = cx + R_light * cos(theta);
            yCircle = cy + R_light * sin(theta);

            if ismember(running, activeIdx)
                plot(xCircle, yCircle, 'r-', 'LineWidth', 1.5);
            elseif ~isempty(inactiveIdx) && ismember(running, inactiveIdx)
                plot(xCircle, yCircle, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
            else
                plot(xCircle, yCircle, 'w-', 'LineWidth', 0.8);
            end
        end
    end
end