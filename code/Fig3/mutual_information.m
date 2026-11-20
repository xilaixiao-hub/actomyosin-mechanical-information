%% ============================================================
% FRAME-RESOLVED SPATIAL MUTUAL INFORMATION
% USING CUMULATIVE TEMPORAL STD MAPS
%
% Required input:
%   TwoSpot_CumulativeTemporalSTD_ProcessedData.mat
%
% Required STD definition:
%   STD_TWO_ALL{e}(:,:,t) = std(rawIntensity(:,:,1:t),0,3)
%
% Primary output:
%   spatialMI{experimentIndex,frameIndex}
%   = bias-corrected spatial MI in nats
%% ============================================================

clearvars;
close all;
clc;

%% ============================================================
% 1. USER SETTINGS
%% ============================================================

inputMAT = fullfile( ...
    'PATH_TO_DATA', ...
    'TwoSpot_CumulativeTemporalSTD_NewExperiments', ...
    'TwoSpot_CumulativeTemporalSTD_ProcessedData.mat');

outputDir = fullfile( ...
    'PATH_TO_DATA', ...
    'TwoSpot_SpatialMI_CumulativeSTD_NewExperiments');



% inputMAT = fullfile( ...
%     'PATH_TO_DATA', ...
%     'TwoSpot_CumulativeTemporalSTD', ...
%     'TwoSpot_CumulativeTemporalSTD_ProcessedData.mat');
% 
% outputDir = fullfile( ...
%     'PATH_TO_DATA', ...
%     'TwoSpot_SpatialMI_CumulativeSTD');


if ~isfolder(outputDir)
    mkdir(outputDir);
end

framesToAnalyze = [];

roiRadiusMultiplier = 2;
maximumAllowedRadiusDifferenceFraction = 0.10;

localGridSpacing_um = [];
maximumLocalGridPointsPerDimension = 301;

binningMode = 'quantile';      % 'quantile' or 'uniform'
numberBins = 6;
minimumPairedPixels = 100;

sampleNormalization = 'none';  % 'none', 'jointZScore', 'separateZScore'

spatialSubsample = 2;

numberSpatialPermutations = 500;
minimumShiftPixels = 5;
rng(1);

frameInterval_s = 10;

saveExampleJointDistributions = true;
exampleExperimentIndex = 1;
exampleFrame = 10;

FZ = 18;
LZ = 2.5;

rgbPurple = [135,85,155]/255;
rgbLightPurple = [205,184,214]/255;
rgbShadePurple = [225,214,231]/255;
rgbGray = [110,110,110]/255;
rgbBlack = [30,30,30]/255;

figurePosition = [100,100,850,680];

saveMAT = true;
saveExcel = true;
saveFigures = true;

%% ============================================================
% 2. LOAD AND VERIFY CUMULATIVE STD DATA
%% ============================================================

if ~isfile(inputMAT)
    error('Input MAT file does not exist:\n%s',inputMAT);
end

S = load(inputMAT);

% requiredVariables = { ...
%     'STD_TWO_ALL','Grid_TWO_ALL','CENTER1','CENTER2', ...
%     'LRadius1','LRadius2','Experiments_Two'};
% 
% for k = 1:numel(requiredVariables)
%     if ~isfield(S,requiredVariables{k})
%         error('Required variable "%s" is missing.',requiredVariables{k});
%     end
% end

requiredVariables = { ...
    'STD_TWO_ALL', ...
    'Grid_TWO_ALL', ...
    'Experiments_Two'};

for k = 1:numel(requiredVariables)

    if ~isfield(S,requiredVariables{k})

        error( ...
            'Required variable "%s" is missing.', ...
            requiredVariables{k});
    end
end



if ~isfield(S,'stdMapMode')
    error(['The MAT file has no stdMapMode field.\n', ...
        'Use the cumulative temporal STD MAT file.']);
end

loadedSTDMode = string(S.stdMapMode);

if ~strcmpi(loadedSTDMode,"cumulativeTemporal")
    error(['Wrong STD dataset selected.\n', ...
        'Expected cumulativeTemporal, but found: %s'],loadedSTDMode);
end

fprintf('\nLoaded STD mode: %s\n',loadedSTDMode);

if isfield(S,'temporalWindowDefinition')
    fprintf('STD definition: %s\n',string(S.temporalWindowDefinition));
end

STD_TWO_ALL = S.STD_TWO_ALL;
Grid_TWO_ALL = S.Grid_TWO_ALL;

% CENTER1 = double(S.CENTER1);
% CENTER2 = double(S.CENTER2);
% 
% LRadius1 = double(S.LRadius1(:));
% LRadius2 = double(S.LRadius2(:));


%% ---------------- Manual analysis-region geometry ----------------
%
% One row per experiment, in the same order as Experiments_Two.
% Units must be micrometers.

CENTER1 = [ ...
    160, 300; ...
    160, 300; ...
    160, 300];

CENTER2 = [ ...
    160, 100; ...
    160, 100; ...
    160, 100];

LRadius1 = [ ...
    20.1; ...
    20.1; ...
    20.1];

LRadius2 = [ ...
    20.1; ...
    20.1; ...
    20.1];

Experiments_Two = double(S.Experiments_Two(:));
numberExperiments = numel(Experiments_Two);

if ~iscell(STD_TWO_ALL)
    error('STD_TWO_ALL must be a cell array.');
end

if numel(STD_TWO_ALL) ~= numberExperiments
    error('STD_TWO_ALL and Experiments_Two have different lengths.');
end

fprintf('Loaded experiments: %s\n',num2str(Experiments_Two'));

%% ============================================================
% 3. DETERMINE FRAME RANGE
%% ============================================================

numberFramesPerExperiment = zeros(numberExperiments,1);

for e = 1:numberExperiments
    currentStack = squeeze(STD_TWO_ALL{e});

    if isempty(currentStack)
        numberFramesPerExperiment(e) = 0;
    elseif ismatrix(currentStack)
        numberFramesPerExperiment(e) = 1;
    elseif ndims(currentStack) == 3
        numberFramesPerExperiment(e) = size(currentStack,3);
    else
        error('Unsupported stack size for experiment %d.',Experiments_Two(e));
    end
end

maximumFrame = max(numberFramesPerExperiment);

if maximumFrame < 1
    error('No valid cumulative STD maps were found.');
end

if isempty(framesToAnalyze)
    framesToAnalyze = 1:maximumFrame;
else
    framesToAnalyze = unique(round(framesToAnalyze(:)'));
    framesToAnalyze = framesToAnalyze( ...
        framesToAnalyze >= 1 & framesToAnalyze <= maximumFrame);
end

if isempty(framesToAnalyze)
    error('No valid frames remain in framesToAnalyze.');
end

numberAnalyzedFrames = numel(framesToAnalyze);

%% ============================================================
% 4. PREALLOCATE RESULTS
%% ============================================================

MI_Raw_bits = nan(numberExperiments,maximumFrame);
MI_NullMean_bits = nan(numberExperiments,maximumFrame);
MI_Excess_bits = nan(numberExperiments,maximumFrame);
MI_BiasCorrected_bits = nan(numberExperiments,maximumFrame);
MI_PermutationP = nan(numberExperiments,maximumFrame);

NumberPairedPixels = nan(numberExperiments,maximumFrame);
PearsonR = nan(numberExperiments,maximumFrame);

LocalGridSpacingUsed_um = nan(numberExperiments,1);
ROIRadiusUsed_um = nan(numberExperiments,1);
LocalGridSize = nan(numberExperiments,2);

PermutationMI_ALL = cell(numberExperiments,maximumFrame);
JointProbability_ALL = cell(numberExperiments,maximumFrame);
IndependentProbability_ALL = cell(numberExperiments,maximumFrame);
EdgesX_ALL = cell(numberExperiments,maximumFrame);
EdgesY_ALL = cell(numberExperiments,maximumFrame);

ConditionPassed = false(numberExperiments,1);
ConditionMessage = strings(numberExperiments,1);

%% ============================================================
% 5. FRAME-RESOLVED SPATIAL MI
%% ============================================================

fprintf('\n============================================\n');
fprintf('CUMULATIVE-STD SPATIAL MUTUAL INFORMATION\n');
fprintf('============================================\n');

for e = 1:numberExperiments

    experimentNumber = Experiments_Two(e);
    fprintf('\nExperiment %d\n',experimentNumber);

    currentSTDStack = squeeze(STD_TWO_ALL{e});
    currentGrid = Grid_TWO_ALL{e};

    if isempty(currentSTDStack) || isempty(currentGrid)
        ConditionMessage(e) = "Missing STD stack or grid";
        continue;
    end

    try
        [XGrid,YGrid] = extractGrid(currentGrid);

        fprintf(['Experiment %d spatial range:\n', ...
    '  x = %.3f to %.3f\n', ...
    '  y = %.3f to %.3f\n'], ...
    experimentNumber, ...
    min(XGrid(:)),max(XGrid(:)), ...
    min(YGrid(:)),max(YGrid(:)));


    catch ME
        ConditionMessage(e) = "Invalid spatial grid";
        warning('Experiment %d: %s',experimentNumber,ME.message);
        continue;
    end

    if ismatrix(currentSTDStack)
        currentSTDStack = reshape( ...
            currentSTDStack,size(currentSTDStack,1),size(currentSTDStack,2),1);
    end

    if size(currentSTDStack,1) == size(XGrid,2) && ...
            size(currentSTDStack,2) == size(XGrid,1) && ...
            ~(size(currentSTDStack,1) == size(XGrid,1) && ...
              size(currentSTDStack,2) == size(XGrid,2))

        currentSTDStack = permute(currentSTDStack,[2,1,3]);
    end

    if size(currentSTDStack,1) ~= size(XGrid,1) || ...
            size(currentSTDStack,2) ~= size(XGrid,2)

        ConditionMessage(e) = "STD-map/grid size mismatch";
        warning('Experiment %d has a stack/grid mismatch.',experimentNumber);
        continue;
    end

    center1 = CENTER1(e,:);
    center2 = CENTER2(e,:);

    R1 = LRadius1(e);
    R2 = LRadius2(e);

    if any(~isfinite([center1 center2 R1 R2])) || R1 <= 0 || R2 <= 0
        ConditionMessage(e) = "Invalid center or radius";
        continue;
    end

    meanRadius = mean([R1 R2]);
    radiusDifferenceFraction = abs(R1-R2)/meanRadius;

    if radiusDifferenceFraction > maximumAllowedRadiusDifferenceFraction
        warning(['Experiment %d has different radii: ', ...
            'R1 = %.3f um, R2 = %.3f um. Using their mean.'], ...
            experimentNumber,R1,R2);
    end

    roiRadius = roiRadiusMultiplier*meanRadius;
    ROIRadiusUsed_um(e) = roiRadius;

    dxValues = abs(diff(XGrid(1,:)));
    dyValues = abs(diff(YGrid(:,1)));

    dxValues = dxValues(isfinite(dxValues) & dxValues > 0);
    dyValues = dyValues(isfinite(dyValues) & dyValues > 0);

    if isempty(dxValues) || isempty(dyValues)
        ConditionMessage(e) = "Invalid grid spacing";
        continue;
    end

    nativeSpacing = max(median(dxValues),median(dyValues));

    if isempty(localGridSpacing_um)
        localSpacing = nativeSpacing;
    else
        localSpacing = localGridSpacing_um;
    end

    estimatedPoints = floor(2*roiRadius/localSpacing)+1;

    if estimatedPoints > maximumLocalGridPointsPerDimension
        localSpacing = 2*roiRadius/(maximumLocalGridPointsPerDimension-1);
    end

    LocalGridSpacingUsed_um(e) = localSpacing;

    relativeCoordinate = -roiRadius:localSpacing:roiRadius;
    [relativeX,relativeY] = meshgrid(relativeCoordinate,relativeCoordinate);

    circularMask = relativeX.^2 + relativeY.^2 <= roiRadius^2;

    subsampleMask = false(size(circularMask));
    subsampleMask(1:spatialSubsample:end,1:spatialSubsample:end) = true;

    analysisMask = circularMask & subsampleMask;

    LocalGridSize(e,:) = size(relativeX);

    queryX1 = center1(1)+relativeX;
    queryY1 = center1(2)+relativeY;

    queryX2 = center2(1)+relativeX;
    queryY2 = center2(2)+relativeY;

    validFrames = framesToAnalyze( ...
        framesToAnalyze <= size(currentSTDStack,3));

    if isempty(validFrames)
        ConditionMessage(e) = "No selected frames available";
        continue;
    end

    ConditionPassed(e) = true;
    ConditionMessage(e) = "Passed";

    for f = validFrames

        currentSTDMap = double(currentSTDStack(:,:,f));

        % Frame 1 of cumulative temporal STD is identically zero.
        if f == 1 || std(currentSTDMap(:),'omitnan') <= eps

            MI_Raw_bits(e,f) = 0;
            MI_NullMean_bits(e,f) = 0;
            MI_Excess_bits(e,f) = 0;
            MI_BiasCorrected_bits(e,f) = 0;
            MI_PermutationP(e,f) = 1;

            fprintf('  Frame %2d: constant cumulative STD; MI set to 0.\n',f);
            continue;
        end

        localMap1 = interp2( ...
            XGrid,YGrid,currentSTDMap,queryX1,queryY1,'linear',NaN);

        localMap2 = interp2( ...
            XGrid,YGrid,currentSTDMap,queryX2,queryY2,'linear',NaN);

        validPairMask = ...
            analysisMask & isfinite(localMap1) & isfinite(localMap2);

        X_samples = localMap1(validPairMask);
        Y_samples = localMap2(validPairMask);

        NumberPairedPixels(e,f) = numel(X_samples);

        if numel(X_samples) < minimumPairedPixels
            warning('Experiment %d frame %d has too few paired pixels.', ...
                experimentNumber,f);
            continue;
        end

        [X_samples,Y_samples] = normalizePairedSamples( ...
            X_samples,Y_samples,sampleNormalization);

        if isempty(X_samples)
            continue;
        end

        [observedMI,Pxy,Pindependent,edgesX,edgesY] = ...
            histogramMutualInformationFixedBins( ...
                X_samples,Y_samples,numberBins,binningMode);

        MI_Raw_bits(e,f) = observedMI;

        JointProbability_ALL{e,f} = Pxy;
        IndependentProbability_ALL{e,f} = Pindependent;
        EdgesX_ALL{e,f} = edgesX;
        EdgesY_ALL{e,f} = edgesY;

        if std(X_samples) > eps && std(Y_samples) > eps
            PearsonR(e,f) = corr( ...
                X_samples,Y_samples,'Type','Pearson','Rows','complete');
        end

        permutationValues = nan(numberSpatialPermutations,1);

        maxRowShift = floor(size(localMap2,1)/2);
        maxColumnShift = floor(size(localMap2,2)/2);

        for p = 1:numberSpatialPermutations

            [rowShift,columnShift] = generateValidSpatialShift( ...
                maxRowShift,maxColumnShift,minimumShiftPixels);

            shiftedMap2 = circshift(localMap2,[rowShift,columnShift]);

            shiftedMask = ...
                analysisMask & isfinite(localMap1) & isfinite(shiftedMap2);

            shiftedX = localMap1(shiftedMask);
            shiftedY = shiftedMap2(shiftedMask);

            if numel(shiftedX) < minimumPairedPixels
                continue;
            end

            [shiftedX,shiftedY] = normalizePairedSamples( ...
                shiftedX,shiftedY,sampleNormalization);

            if isempty(shiftedX)
                continue;
            end

            permutationValues(p) = histogramMutualInformationFixedBins( ...
                shiftedX,shiftedY,numberBins,binningMode,edgesX,edgesY);
        end

        permutationValues = permutationValues(isfinite(permutationValues));
        PermutationMI_ALL{e,f} = permutationValues;

        if isempty(permutationValues)
            continue;
        end

        nullMean = mean(permutationValues,'omitnan');
        excessMI = observedMI-nullMean;
        correctedMI = max(excessMI,0);

        permutationP = ...
            (1+sum(permutationValues >= observedMI))/ ...
            (numel(permutationValues)+1);

        MI_NullMean_bits(e,f) = nullMean;
        MI_Excess_bits(e,f) = excessMI;
        MI_BiasCorrected_bits(e,f) = correctedMI;
        MI_PermutationP(e,f) = permutationP;

        fprintf(['  Frame %2d: raw = %.4f bits, null = %.4f, ', ...
            'corrected = %.4f, p = %.4g, N = %d\n'], ...
            f,observedMI,nullMean,correctedMI,permutationP,numel(X_samples));
    end
end

%% ============================================================
% 6. CONVERT TO NATS AND CELL ARRAYS
%% ============================================================

MI_Raw_nats = MI_Raw_bits*log(2);
MI_NullMean_nats = MI_NullMean_bits*log(2);
MI_Excess_nats = MI_Excess_bits*log(2);
MI_BiasCorrected_nats = MI_BiasCorrected_bits*log(2);

spatialMI = numericMatrixToCell(MI_BiasCorrected_nats);
spatialMI_bits = numericMatrixToCell(MI_BiasCorrected_bits);
spatialMI_raw_bits = numericMatrixToCell(MI_Raw_bits);
spatialMI_null_bits = numericMatrixToCell(MI_NullMean_bits);
spatialMI_excess_bits = numericMatrixToCell(MI_Excess_bits);
spatialMI_p = numericMatrixToCell(MI_PermutationP);

%% ============================================================
% 7. GROUP STATISTICS
%% ============================================================

MeanMI_Raw_bits = nan(1,maximumFrame);
SEM_MI_Raw_bits = nan(1,maximumFrame);

MeanMI_Corrected_bits = nan(1,maximumFrame);
SD_MI_Corrected_bits = nan(1,maximumFrame);
SEM_MI_Corrected_bits = nan(1,maximumFrame);

MeanMI_Corrected_nats = nan(1,maximumFrame);
SD_MI_Corrected_nats = nan(1,maximumFrame);
SEM_MI_Corrected_nats = nan(1,maximumFrame);

NumberValidExperimentsByFrame = zeros(1,maximumFrame);

for f = framesToAnalyze

    rawValues = MI_Raw_bits(:,f);
    correctedBits = MI_BiasCorrected_bits(:,f);
    correctedNats = MI_BiasCorrected_nats(:,f);

    validRaw = isfinite(rawValues);
    validCorrected = isfinite(correctedBits);

    NumberValidExperimentsByFrame(f) = nnz(validCorrected);

    if any(validRaw)
        MeanMI_Raw_bits(f) = mean(rawValues(validRaw));
        SEM_MI_Raw_bits(f) = std(rawValues(validRaw),0)/sqrt(nnz(validRaw));
    end

    if any(validCorrected)
        MeanMI_Corrected_bits(f) = mean(correctedBits(validCorrected));
        SD_MI_Corrected_bits(f) = std(correctedBits(validCorrected),0);
        SEM_MI_Corrected_bits(f) = ...
            SD_MI_Corrected_bits(f)/sqrt(nnz(validCorrected));

        MeanMI_Corrected_nats(f) = mean(correctedNats(validCorrected));
        SD_MI_Corrected_nats(f) = std(correctedNats(validCorrected),0);
        SEM_MI_Corrected_nats(f) = ...
            SD_MI_Corrected_nats(f)/sqrt(nnz(validCorrected));
    end
end

meanNullByFrame_bits = mean(MI_NullMean_bits,1,'omitnan');

%% ============================================================
% 8. TIME AXIS
%% ============================================================

if isempty(frameInterval_s)
    xVector = framesToAnalyze;
    xAxisLabel = 'Frame';
else
    xVector = (framesToAnalyze-1)*frameInterval_s;
    xAxisLabel = 'Time (s)';
end

%% ============================================================
% 9. FIGURES
%% ============================================================

plotData = MI_BiasCorrected_nats(:,framesToAnalyze);
meanMI = MeanMI_Corrected_nats(framesToAnalyze);
stdMI = SD_MI_Corrected_nats(framesToAnalyze);

lowerSD = max(meanMI-stdMI,0);
upperSD = meanMI+stdMI;

figIndividual = figure('Color','w','Position',figurePosition);
ax = axes(figIndividual);
hold(ax,'on');

validShade = isfinite(lowerSD) & isfinite(upperSD);

if nnz(validShade) >= 2
    fill(ax, ...
        [xVector(validShade),fliplr(xVector(validShade))], ...
        [upperSD(validShade),fliplr(lowerSD(validShade))], ...
        rgbShadePurple,'FaceAlpha',0.60,'EdgeColor','none', ...
        'DisplayName','Mean \pm SD');
end

for e = 1:numberExperiments
    plot(ax,xVector,plotData(e,:),'-o', ...
        'Color',rgbLightPurple,'LineWidth',1.4, ...
        'MarkerSize',4.5,'MarkerFaceColor',rgbLightPurple, ...
        'MarkerEdgeColor','none','HandleVisibility','off');
end

plot(ax,xVector,meanMI,'-o', ...
    'Color',rgbPurple,'LineWidth',3.5, ...
    'MarkerSize',6.5,'MarkerFaceColor',rgbPurple, ...
    'MarkerEdgeColor','none','DisplayName','Mean');

xlabel(ax,xAxisLabel);
ylabel(ax,'Spatial mutual information (nats)');
title(ax,'Cumulative temporal-STD spatial MI','FontWeight','normal');
formatPaperAxes(ax,FZ,rgbBlack);
legend(ax,'Location','northwest','Box','off');

if numel(xVector) >= 2
    xlim(ax,[min(xVector),max(xVector)]);
end

yMax = max([upperSD(isfinite(upperSD)),plotData(isfinite(plotData))']);

if isempty(yMax) || ~isfinite(yMax) || yMax <= 0
    yMax = 1;
end

ylim(ax,[0,1.10*yMax]);
hold(ax,'off');

figAverage = figure('Color','w','Position',figurePosition);
ax = axes(figAverage);
hold(ax,'on');

semValues = SEM_MI_Corrected_nats(framesToAnalyze);
validSEM = isfinite(meanMI) & isfinite(semValues);

if nnz(validSEM) >= 2

    lowerSEM = max(meanMI(validSEM)-semValues(validSEM),0);
    upperSEM = meanMI(validSEM)+semValues(validSEM);

    fill(ax, ...
        [xVector(validSEM),fliplr(xVector(validSEM))], ...
        [lowerSEM,fliplr(upperSEM)], ...
        rgbShadePurple,'FaceAlpha',0.65,'EdgeColor','none');
end

plot(ax,xVector,meanMI,'-o', ...
    'Color',rgbPurple,'LineWidth',3, ...
    'MarkerSize',7,'MarkerFaceColor',rgbPurple, ...
    'MarkerEdgeColor','none');

xlabel(ax,xAxisLabel);
ylabel(ax,'Spatial mutual information (nats)');
title(ax,'Mean bias-corrected cumulative-STD MI','FontWeight','normal');
formatPaperAxes(ax,FZ,rgbBlack);

if numel(xVector) >= 2
    xlim(ax,[min(xVector),max(xVector)]);
end

hold(ax,'off');

figRawNull = figure('Color','w','Position',figurePosition);
ax = axes(figRawNull);
hold(ax,'on');

plot(ax,xVector,MeanMI_Raw_bits(framesToAnalyze)*log(2),'-', ...
    'Color',rgbPurple,'LineWidth',2.8,'DisplayName','Observed MI');

plot(ax,xVector,meanNullByFrame_bits(framesToAnalyze)*log(2),'--', ...
    'Color',rgbGray,'LineWidth',2.2,'DisplayName','Spatial-shift null');

xlabel(ax,xAxisLabel);
ylabel(ax,'Mutual information (nats)');
title(ax,'Observed and null cumulative-STD MI','FontWeight','normal');
legend(ax,'Location','best','Box','off');
formatPaperAxes(ax,FZ,rgbBlack);

if numel(xVector) >= 2
    xlim(ax,[min(xVector),max(xVector)]);
end

hold(ax,'off');

figHeatmap = figure('Color','w','Position',[100,100,900,560]);
ax = axes(figHeatmap);

imagesc(ax,xVector,1:numberExperiments, ...
    MI_BiasCorrected_nats(:,framesToAnalyze));

set(ax,'YDir','normal');

colormap(ax,purpleProbabilityMap(256,rgbPurple));

cb = colorbar(ax);
cb.Label.String = 'Bias-corrected spatial MI (nats)';

xlabel(ax,xAxisLabel);
ylabel(ax,'Two-spot experiment');

yticks(ax,1:numberExperiments);
yticklabels(ax,string(Experiments_Two));

title(ax,'Cumulative temporal-STD spatial MI','FontWeight','normal');
formatPaperAxes(ax,FZ,rgbBlack);

%% ============================================================
% 10. EXAMPLE DIAGNOSTIC FIGURES
%% ============================================================

figExampleMaps = [];
figExampleJoint = [];
figExampleNull = [];

if saveExampleJointDistributions && ...
        exampleExperimentIndex >= 1 && ...
        exampleExperimentIndex <= numberExperiments && ...
        exampleFrame >= 1 && ...
        exampleFrame <= maximumFrame && ...
        isfinite(MI_Raw_bits(exampleExperimentIndex,exampleFrame))

    e = exampleExperimentIndex;

    currentSTDStack = squeeze(STD_TWO_ALL{e});

    if ismatrix(currentSTDStack)
        currentSTDStack = reshape( ...
            currentSTDStack,size(currentSTDStack,1),size(currentSTDStack,2),1);
    end

    [XGrid,YGrid] = extractGrid(Grid_TWO_ALL{e});

    if size(currentSTDStack,1) == size(XGrid,2) && ...
            size(currentSTDStack,2) == size(XGrid,1) && ...
            ~(size(currentSTDStack,1) == size(XGrid,1) && ...
              size(currentSTDStack,2) == size(XGrid,2))

        currentSTDStack = permute(currentSTDStack,[2,1,3]);
    end

    radius = ROIRadiusUsed_um(e);
    spacing = LocalGridSpacingUsed_um(e);

    relativeCoordinate = -radius:spacing:radius;
    [relativeX,relativeY] = meshgrid(relativeCoordinate,relativeCoordinate);

    circularMask = relativeX.^2 + relativeY.^2 <= radius^2;

    currentMap = double(currentSTDStack(:,:,exampleFrame));

    map1 = interp2(XGrid,YGrid,currentMap, ...
        CENTER1(e,1)+relativeX,CENTER1(e,2)+relativeY,'linear',NaN);

    map2 = interp2(XGrid,YGrid,currentMap, ...
        CENTER2(e,1)+relativeX,CENTER2(e,2)+relativeY,'linear',NaN);

    map1(~circularMask) = NaN;
    map2(~circularMask) = NaN;

    figExampleMaps = figure('Color','w','Position',[100,100,1050,480]);
    tiledlayout(figExampleMaps,1,2,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile;
    imagesc(ax1,relativeCoordinate,relativeCoordinate,map1);
    axis(ax1,'image');
    set(ax1,'YDir','normal');
    colormap(ax1,purpleProbabilityMap(256,rgbPurple));
    colorbar(ax1);
    xlabel(ax1,'\Deltax (\mum)');
    ylabel(ax1,'\Deltay (\mum)');
    title(ax1,'Spot 1 cumulative STD map','FontWeight','normal');
    formatPaperAxes(ax1,FZ,rgbBlack);

    ax2 = nexttile;
    imagesc(ax2,relativeCoordinate,relativeCoordinate,map2);
    axis(ax2,'image');
    set(ax2,'YDir','normal');
    colormap(ax2,purpleProbabilityMap(256,rgbPurple));
    colorbar(ax2);
    xlabel(ax2,'\Deltax (\mum)');
    ylabel(ax2,'\Deltay (\mum)');
    title(ax2,'Spot 2 cumulative STD map','FontWeight','normal');
    formatPaperAxes(ax2,FZ,rgbBlack);

    sgtitle(figExampleMaps,sprintf( ...
        'Experiment %d, frame %d, corrected MI = %.4f nats', ...
        Experiments_Two(e),exampleFrame,MI_BiasCorrected_nats(e,exampleFrame)), ...
        'FontWeight','normal');

    Pxy = JointProbability_ALL{e,exampleFrame};
    Pindependent = IndependentProbability_ALL{e,exampleFrame};
    edgesX = EdgesX_ALL{e,exampleFrame};
    edgesY = EdgesY_ALL{e,exampleFrame};

    if ~isempty(Pxy)

        xCenters = 0.5*(edgesX(1:end-1)+edgesX(2:end));
        yCenters = 0.5*(edgesY(1:end-1)+edgesY(2:end));

        figExampleJoint = figure('Color','w','Position',[100,100,1040,480]);
        tiledlayout(figExampleJoint,1,2,'TileSpacing','compact','Padding','compact');

        ax1 = nexttile;
        imagesc(ax1,xCenters,yCenters,Pxy');
        axis(ax1,'xy');
        axis(ax1,'square');
        colormap(ax1,purpleProbabilityMap(256,rgbPurple));
        colorbar(ax1);
        xlabel(ax1,'Spot 1 cumulative STD');
        ylabel(ax1,'Spot 2 cumulative STD');
        title(ax1,'Observed joint probability','FontWeight','normal');
        formatPaperAxes(ax1,FZ,rgbBlack);

        ax2 = nexttile;
        imagesc(ax2,xCenters,yCenters,Pindependent');
        axis(ax2,'xy');
        axis(ax2,'square');
        colormap(ax2,purpleProbabilityMap(256,rgbPurple));
        colorbar(ax2);
        xlabel(ax2,'Spot 1 cumulative STD');
        ylabel(ax2,'Spot 2 cumulative STD');
        title(ax2,'Independent expectation','FontWeight','normal');
        formatPaperAxes(ax2,FZ,rgbBlack);
    end

    nullValues = PermutationMI_ALL{e,exampleFrame};

    if ~isempty(nullValues)

        figExampleNull = figure('Color','w','Position',figurePosition);
        ax = axes(figExampleNull);
        hold(ax,'on');

        histogram(ax,nullValues,30, ...
            'FaceColor',rgbLightPurple,'EdgeColor','none');

        xline(ax,MI_Raw_bits(e,exampleFrame),'-', ...
            sprintf('Observed = %.3f bits',MI_Raw_bits(e,exampleFrame)), ...
            'Color',rgbPurple,'LineWidth',LZ);

        xlabel(ax,'Spatial-shift MI (bits)');
        ylabel(ax,'Count');

        title(ax,sprintf('Experiment %d, frame %d, p = %.4g', ...
            Experiments_Two(e),exampleFrame,MI_PermutationP(e,exampleFrame)), ...
            'FontWeight','normal');

        formatPaperAxes(ax,FZ,rgbBlack);
        hold(ax,'off');
    end
end

%% ============================================================
% 11. OUTPUT TABLES
%% ============================================================

ExperimentColumn = [];
FrameColumn = [];
TimeColumn_s = [];

RawMI_bits_Column = [];
NullMean_bits_Column = [];
ExcessMI_bits_Column = [];
CorrectedMI_bits_Column = [];
CorrectedMI_nats_Column = [];

PermutationPColumn = [];
PearsonRColumn = [];
NumberPairsColumn = [];

for e = 1:numberExperiments

    for f = framesToAnalyze

        if ~isfinite(MI_Raw_bits(e,f))
            continue;
        end

        ExperimentColumn(end+1,1) = Experiments_Two(e); %#ok<SAGROW>
        FrameColumn(end+1,1) = f; %#ok<SAGROW>

        if isempty(frameInterval_s)
            TimeColumn_s(end+1,1) = NaN; %#ok<SAGROW>
        else
            TimeColumn_s(end+1,1) = (f-1)*frameInterval_s; %#ok<SAGROW>
        end

        RawMI_bits_Column(end+1,1) = MI_Raw_bits(e,f); %#ok<SAGROW>
        NullMean_bits_Column(end+1,1) = MI_NullMean_bits(e,f); %#ok<SAGROW>
        ExcessMI_bits_Column(end+1,1) = MI_Excess_bits(e,f); %#ok<SAGROW>
        CorrectedMI_bits_Column(end+1,1) = MI_BiasCorrected_bits(e,f); %#ok<SAGROW>
        CorrectedMI_nats_Column(end+1,1) = MI_BiasCorrected_nats(e,f); %#ok<SAGROW>
        PermutationPColumn(end+1,1) = MI_PermutationP(e,f); %#ok<SAGROW>
        PearsonRColumn(end+1,1) = PearsonR(e,f); %#ok<SAGROW>
        NumberPairsColumn(end+1,1) = NumberPairedPixels(e,f); %#ok<SAGROW>
    end
end

FrameResolvedTable = table( ...
    ExperimentColumn,FrameColumn,TimeColumn_s, ...
    RawMI_bits_Column,NullMean_bits_Column,ExcessMI_bits_Column, ...
    CorrectedMI_bits_Column,CorrectedMI_nats_Column, ...
    PermutationPColumn,PearsonRColumn,NumberPairsColumn, ...
    'VariableNames',{ ...
    'Experiment','Frame','Time_s','RawSpatialMI_bits', ...
    'SpatialShiftNullMean_bits','ExcessSpatialMI_bits', ...
    'BiasCorrectedSpatialMI_bits','BiasCorrectedSpatialMI_nats', ...
    'PermutationP','PearsonR','NumberPairedSpatialSamples'});

FrameSummaryTable = table( ...
    framesToAnalyze(:),reshape(xVector,[],1), ...
    MeanMI_Raw_bits(framesToAnalyze)', ...
    SEM_MI_Raw_bits(framesToAnalyze)', ...
    meanNullByFrame_bits(framesToAnalyze)', ...
    MeanMI_Corrected_bits(framesToAnalyze)', ...
    SD_MI_Corrected_bits(framesToAnalyze)', ...
    SEM_MI_Corrected_bits(framesToAnalyze)', ...
    MeanMI_Corrected_nats(framesToAnalyze)', ...
    SD_MI_Corrected_nats(framesToAnalyze)', ...
    SEM_MI_Corrected_nats(framesToAnalyze)', ...
    NumberValidExperimentsByFrame(framesToAnalyze)', ...
    'VariableNames',{ ...
    'Frame','TimeOrFrameCoordinate','MeanRawMI_bits','SEMRawMI_bits', ...
    'MeanSpatialShiftNull_bits','MeanCorrectedMI_bits', ...
    'SDCorrectedMI_bits','SEMCorrectedMI_bits', ...
    'MeanCorrectedMI_nats','SDCorrectedMI_nats', ...
    'SEMCorrectedMI_nats','NumberValidExperiments'});

ExperimentSettingsTable = table( ...
    Experiments_Two,numberFramesPerExperiment,ROIRadiusUsed_um, ...
    LocalGridSpacingUsed_um,LocalGridSize(:,1),LocalGridSize(:,2), ...
    ConditionPassed,ConditionMessage, ...
    'VariableNames',{ ...
    'Experiment','NumberAvailableFrames','ROIRadiusUsed_um', ...
    'LocalGridSpacingUsed_um','LocalGridRows','LocalGridColumns', ...
    'ConditionPassed','ConditionMessage'});

%% ============================================================
% 12. SAVE OUTPUTS
%% ============================================================

simpleMIFile = fullfile(outputDir,'spatialMI_CumulativeSTD.mat');

save(simpleMIFile, ...
    'spatialMI','spatialMI_bits','spatialMI_raw_bits', ...
    'spatialMI_null_bits','spatialMI_excess_bits','spatialMI_p', ...
    'MI_BiasCorrected_nats','MI_BiasCorrected_bits', ...
    'MI_Raw_bits','MI_NullMean_bits','MI_Excess_bits', ...
    'MI_PermutationP','Experiments_Two','framesToAnalyze', ...
    'frameInterval_s','loadedSTDMode');

if saveExcel

    outputExcel = fullfile(outputDir,'TwoSpot_SpatialMI_CumulativeSTD.xlsx');

    if isfile(outputExcel)
        delete(outputExcel);
    end

    writetable(FrameResolvedTable,outputExcel,'Sheet','FrameResolvedMI');
    writetable(FrameSummaryTable,outputExcel,'Sheet','FrameSummary');
    writetable(ExperimentSettingsTable,outputExcel,'Sheet','ExperimentSettings');

    writematrix(MI_Raw_bits,outputExcel,'Sheet','RawMI_bits');
    writematrix(MI_NullMean_bits,outputExcel,'Sheet','NullMeanMI_bits');
    writematrix(MI_Excess_bits,outputExcel,'Sheet','ExcessMI_bits');
    writematrix(MI_BiasCorrected_bits,outputExcel,'Sheet','CorrectedMI_bits');
    writematrix(MI_BiasCorrected_nats,outputExcel,'Sheet','CorrectedMI_nats');
    writematrix(MI_PermutationP,outputExcel,'Sheet','PermutationP');
end

if saveMAT

    outputMAT = fullfile(outputDir,'TwoSpot_SpatialMI_CumulativeSTD.mat');

    save(outputMAT, ...
        'spatialMI','spatialMI_bits','spatialMI_raw_bits', ...
        'spatialMI_null_bits','spatialMI_excess_bits','spatialMI_p', ...
        'MI_Raw_bits','MI_NullMean_bits','MI_Excess_bits', ...
        'MI_BiasCorrected_bits','MI_Raw_nats','MI_NullMean_nats', ...
        'MI_Excess_nats','MI_BiasCorrected_nats','MI_PermutationP', ...
        'NumberPairedPixels','PearsonR','PermutationMI_ALL', ...
        'JointProbability_ALL','IndependentProbability_ALL', ...
        'EdgesX_ALL','EdgesY_ALL','MeanMI_Raw_bits', ...
        'SEM_MI_Raw_bits','MeanMI_Corrected_bits', ...
        'SD_MI_Corrected_bits','SEM_MI_Corrected_bits', ...
        'MeanMI_Corrected_nats','SD_MI_Corrected_nats', ...
        'SEM_MI_Corrected_nats','NumberValidExperimentsByFrame', ...
        'Experiments_Two','framesToAnalyze','frameInterval_s', ...
        'ROIRadiusUsed_um','LocalGridSpacingUsed_um','LocalGridSize', ...
        'ConditionPassed','ConditionMessage','FrameResolvedTable', ...
        'FrameSummaryTable','ExperimentSettingsTable','loadedSTDMode', ...
        'roiRadiusMultiplier','maximumAllowedRadiusDifferenceFraction', ...
        'numberBins','binningMode','minimumPairedPixels', ...
        'sampleNormalization','spatialSubsample', ...
        'numberSpatialPermutations','minimumShiftPixels','-v7.3');
end

if saveFigures

    saveVectorAndPNG(figIndividual,outputDir, ...
        'CumulativeSTD_SpatialMI_Individual_Mean_SD');

    saveVectorAndPNG(figAverage,outputDir, ...
        'CumulativeSTD_SpatialMI_Mean_SEM');

    saveVectorAndPNG(figRawNull,outputDir, ...
        'CumulativeSTD_SpatialMI_Raw_Vs_Null');

    exportgraphics(figHeatmap, ...
        fullfile(outputDir,'CumulativeSTD_SpatialMI_Heatmap.png'), ...
        'Resolution',600);

    if ~isempty(figExampleMaps)
        exportgraphics(figExampleMaps, ...
            fullfile(outputDir,'CumulativeSTD_SpatialMI_ExampleMaps.png'), ...
            'Resolution',600);
    end

    if ~isempty(figExampleJoint)
        exportgraphics(figExampleJoint, ...
            fullfile(outputDir,'CumulativeSTD_SpatialMI_ExampleJoint.png'), ...
            'Resolution',600);
    end

    if ~isempty(figExampleNull)
        saveVectorAndPNG(figExampleNull,outputDir, ...
            'CumulativeSTD_SpatialMI_ExampleNull');
    end
end

fprintf('\n============================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('============================================\n');
fprintf('STD mode: %s\n',loadedSTDMode);
fprintf('Experiments: %d\n',numberExperiments);
fprintf('Analyzed frames: %d\n',numberAnalyzedFrames);
fprintf('Primary output: spatialMI in nats\n');
fprintf('Output folder:\n%s\n',outputDir);

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function [XGrid,YGrid] = extractGrid(currentGrid)

    if ~iscell(currentGrid) || numel(currentGrid) < 2
        error('Grid entry must contain X and Y coordinates.');
    end

    XInput = double(currentGrid{1});
    YInput = double(currentGrid{2});

    if isvector(XInput) && isvector(YInput)

        [XGrid,YGrid] = meshgrid(XInput(:)',YInput(:));

    elseif isequal(size(XInput),size(YInput))

        XGrid = XInput;
        YGrid = YInput;

    else

        error('X and Y grid dimensions are incompatible.');
    end
end


function cellArray = numericMatrixToCell(matrixValue)

    cellArray = cell(size(matrixValue));

    for r = 1:size(matrixValue,1)
        for c = 1:size(matrixValue,2)

            if isfinite(matrixValue(r,c))
                cellArray{r,c} = matrixValue(r,c);
            else
                cellArray{r,c} = [];
            end
        end
    end
end


function [X,Y] = normalizePairedSamples(X,Y,normalizationMode)

    X = X(:);
    Y = Y(:);

    valid = isfinite(X) & isfinite(Y);
    X = X(valid);
    Y = Y(valid);

    switch lower(normalizationMode)

        case 'none'

        case 'jointzscore'

            combined = [X;Y];
            mu = mean(combined,'omitnan');
            sigma = std(combined,0,'omitnan');

            if sigma <= eps
                X = [];
                Y = [];
                return;
            end

            X = (X-mu)/sigma;
            Y = (Y-mu)/sigma;

        case 'separatezscore'

            sigmaX = std(X,0,'omitnan');
            sigmaY = std(Y,0,'omitnan');

            if sigmaX <= eps || sigmaY <= eps
                X = [];
                Y = [];
                return;
            end

            X = (X-mean(X,'omitnan'))/sigmaX;
            Y = (Y-mean(Y,'omitnan'))/sigmaY;

        otherwise

            error('Unknown sampleNormalization mode.');
    end
end


function [MI_bits,Pxy,Pindependent,edgesX,edgesY] = ...
        histogramMutualInformationFixedBins( ...
            X,Y,numberBins,binningMode,suppliedEdgesX,suppliedEdgesY)

    X = X(:);
    Y = Y(:);

    valid = isfinite(X) & isfinite(Y);
    X = X(valid);
    Y = Y(valid);

    if numel(X) < 2

        MI_bits = NaN;
        Pxy = [];
        Pindependent = [];
        edgesX = [];
        edgesY = [];
        return;
    end

    if nargin >= 6 && ...
            ~isempty(suppliedEdgesX) && ...
            ~isempty(suppliedEdgesY)

        edgesX = suppliedEdgesX;
        edgesY = suppliedEdgesY;

    else

        switch lower(binningMode)

            case 'quantile'

                q = linspace(0,1,numberBins+1);
                edgesX = quantile(X,q);
                edgesY = quantile(Y,q);

            case 'uniform'

                edgesX = linspace(min(X),max(X),numberBins+1);
                edgesY = linspace(min(Y),max(Y),numberBins+1);

            otherwise

                error('Unknown binningMode.');
        end

        edgesX = makeStrictlyIncreasingEdges(edgesX,X);
        edgesY = makeStrictlyIncreasingEdges(edgesY,Y);
    end

    jointCounts = histcounts2(X,Y,edgesX,edgesY);
    totalCount = sum(jointCounts,'all');

    if totalCount <= 0

        MI_bits = NaN;
        Pxy = [];
        Pindependent = [];
        return;
    end

    Pxy = jointCounts/totalCount;

    Px = sum(Pxy,2);
    Py = sum(Pxy,1);

    Pindependent = Px*Py;

    validProbability = Pxy > 0 & Pindependent > 0;

    MI_bits = sum( ...
        Pxy(validProbability).* ...
        log2(Pxy(validProbability)./Pindependent(validProbability)), ...
        'all');
end


function edges = makeStrictlyIncreasingEdges(edges,data)

    edges = edges(:)';
    data = data(isfinite(data));

    dataMin = min(data);
    dataMax = max(data);
    dataRange = dataMax-dataMin;

    if dataRange <= eps

        centerValue = mean(data);

        edges = linspace( ...
            centerValue-0.5,centerValue+0.5,numel(edges));

        return;
    end

    increment = max(eps(max(abs(data))),1e-10*dataRange);

    for k = 2:numel(edges)

        if edges(k) <= edges(k-1)
            edges(k) = edges(k-1)+increment;
        end
    end

    edges(1) = min(edges(1),dataMin-increment);
    edges(end) = max(edges(end),dataMax+increment);
end


function [rowShift,columnShift] = ...
        generateValidSpatialShift( ...
            maximumRowShift,maximumColumnShift,minimumShiftPixels)

    for attempt = 1:1000

        rowShift = randi([-maximumRowShift,maximumRowShift]);
        columnShift = randi([-maximumColumnShift,maximumColumnShift]);

        if hypot(rowShift,columnShift) >= minimumShiftPixels
            return;
        end
    end

    rowShift = min(minimumShiftPixels,maximumRowShift);
    columnShift = min(minimumShiftPixels,maximumColumnShift);
end


function formatPaperAxes(axisHandle,fontSize,blackColor)

    set(axisHandle, ...
        'FontName','Arial', ...
        'FontSize',fontSize, ...
        'LineWidth',1.2, ...
        'TickDir','in', ...
        'TickLength',[0.018,0.018], ...
        'Box','on', ...
        'XColor',blackColor, ...
        'YColor',blackColor, ...
        'Layer','top');

    grid(axisHandle,'off');
end


function cmap = purpleProbabilityMap(numberColors,mainPurple)

    whiteColor = [1,1,1];
    darkPurple = 0.40*mainPurple;

    anchorColors = [ ...
        whiteColor; ...
        0.75*whiteColor+0.25*mainPurple; ...
        mainPurple; ...
        darkPurple];

    anchorPositions = linspace(0,1,size(anchorColors,1));
    requestedPositions = linspace(0,1,numberColors);

    cmap = interp1( ...
        anchorPositions,anchorColors,requestedPositions,'pchip');

    cmap = min(max(cmap,0),1);
end


function saveVectorAndPNG(figureHandle,outputDirectory,baseName)

    if isempty(figureHandle) || ~isgraphics(figureHandle)
        return;
    end

    set(figureHandle,'Renderer','painters');

    exportgraphics(figureHandle, ...
        fullfile(outputDirectory,[baseName '.pdf']), ...
        'ContentType','vector');

    exportgraphics(figureHandle, ...
        fullfile(outputDirectory,[baseName '.eps']), ...
        'ContentType','vector');

    exportgraphics(figureHandle, ...
        fullfile(outputDirectory,[baseName '.png']), ...
        'Resolution',600);
end
