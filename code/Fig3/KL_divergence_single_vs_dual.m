%% ============================================================
% KL WITH MANUAL CONTROL REGIONS
% CONSISTENT WITH THE ORIGINAL PIXEL-VOLATILITY KL SCRIPT
%% ============================================================

clearvars;
close all;
clc;

%% ============================================================
% 1. USER SETTINGS
%% ============================================================

Experiments_Single = [24 25 26 28];
Experiments_Two    = [6 8 9 10 11 13];

singleDataBasePath = fullfile( ...
    'PATH_TO_DATA', ...
    'Experiments', ...
    'SingleLightSpot_Data&PIV');

singleFolderPrefix = 'set';

twoSpotMAT = fullfile( ...
    'PATH_TO_DATA', ...
    'TwoSpot_CumulativeTemporalSTD', ...
    'TwoSpot_CumulativeTemporalSTD_ProcessedData.mat');

outputDir = fullfile( ...
    'PATH_TO_DATA', ...
    'TwoSpot_ActivatedVsControl_KL_OriginalMethod');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

numberFramesRequested = 18;
framesToAnalyze = [];
frameInterval_s = 10;

singleReferenceRadiusMultiplier = 2;
activatedROIRadiusMultiplier = 2;

normalizationLowerPercentile = 1;
normalizationUpperPercentile = 99;

removeOutliersBeforeKL = false;
outlierIQRMultiplier = 2;

evaluationPoints = linspace(-3, 3, 256);
kdeEpsilon = 1e-12;
minimumSamplesPerDistribution = 20;

calculateControl = true;
controlSelectionDisplay = 'temporalMean';
controlSelectionFrame = 10;
reuseSavedControlROIs = true;

controlROIFile = fullfile( ...
    outputDir, ...
    'Selected_Nonactivated_Control_ROIs_KL.mat');

fixedControlSquareSide_um = [];
controlExclusionMultiplier = 2;

saveKLDiagnosticPlots = false;
diagnosticExperiments = [];
diagnosticFrames = [];

FZ = 18;

rgb_Purple      = [135, 85, 155]/255;
rgb_LightPurple = [202, 178, 214]/255;
rgb_ShadePurple = [225, 211, 232]/255;

rgb_Gray      = [105, 105, 105]/255;
rgb_LightGray = [185, 185, 185]/255;
rgb_ShadeGray = [225, 225, 225]/255;

rgb_Black  = [30, 30, 30]/255;
rgb_Cyan   = [55, 185, 205]/255;
rgb_Orange = [235, 145, 40]/255;

figurePosition = [100, 100, 880, 680];

saveMAT = true;
saveExcel = true;
saveFigures = true;

%% ============================================================
% 2. LOAD CUMULATIVE TWO-SPOT STD DATA
%% ============================================================

if ~isfile(twoSpotMAT)
    error('Two-spot MAT file does not exist:\n%s', twoSpotMAT);
end

TwoSpotData = load(twoSpotMAT);

requiredTwoSpotVariables = { ...
    'STD_TWO_ALL', ...
    'Grid_TWO_ALL', ...
    'CENTER1', ...
    'CENTER2', ...
    'LRadius1', ...
    'LRadius2', ...
    'Experiments_Two'};

for variableIndex = 1:numel(requiredTwoSpotVariables)
    variableName = requiredTwoSpotVariables{variableIndex};
    if ~isfield(TwoSpotData, variableName)
        error('Two-spot MAT file is missing "%s".', variableName);
    end
end

if isfield(TwoSpotData, 'stdMapMode')
    loadedMode = string(TwoSpotData.stdMapMode);
    if ~strcmpi(loadedMode, "cumulativeTemporal")
        error('Expected cumulativeTemporal, but found: %s', loadedMode);
    end
else
    warning(['The MAT file has no stdMapMode field. Verify that ', ...
             'STD_TWO_ALL contains cumulative temporal STD maps.']);
end

Saved_STD_TWO_ALL = TwoSpotData.STD_TWO_ALL;
Saved_Grid_TWO_ALL = TwoSpotData.Grid_TWO_ALL;
Saved_CENTER1 = double(TwoSpotData.CENTER1);
Saved_CENTER2 = double(TwoSpotData.CENTER2);
Saved_LRadius1 = double(TwoSpotData.LRadius1(:));
Saved_LRadius2 = double(TwoSpotData.LRadius2(:));
Saved_Experiments_Two = double(TwoSpotData.Experiments_Two(:));

Experiments_Two = double(Experiments_Two(:));
numberTwoSpotExperiments = numel(Experiments_Two);

STD_TWO_ALL = cell(numberTwoSpotExperiments, 1);
Grid_TWO_ALL = cell(numberTwoSpotExperiments, 1);
CENTER1 = nan(numberTwoSpotExperiments, 2);
CENTER2 = nan(numberTwoSpotExperiments, 2);
LRadius1 = nan(numberTwoSpotExperiments, 1);
LRadius2 = nan(numberTwoSpotExperiments, 1);

for experimentIndex = 1:numberTwoSpotExperiments
    experimentNumber = Experiments_Two(experimentIndex);
    savedIndex = find(Saved_Experiments_Two == experimentNumber, 1, 'first');

    if isempty(savedIndex)
        error('Experiment %d is not present in the saved MAT file.', experimentNumber);
    end

    currentStack = squeeze(Saved_STD_TWO_ALL{savedIndex});
    if ismatrix(currentStack)
        currentStack = reshape(currentStack, size(currentStack,1), size(currentStack,2), 1);
    end

    [XGrid, YGrid, xVector, yVector] = ...
        extractSpatialGrid(Saved_Grid_TWO_ALL{savedIndex});

    if size(currentStack,1) == size(XGrid,2) && ...
            size(currentStack,2) == size(XGrid,1) && ...
            ~(size(currentStack,1) == size(XGrid,1) && ...
              size(currentStack,2) == size(XGrid,2))
        currentStack = permute(currentStack, [2 1 3]);
    end

    if size(currentStack,1) ~= size(XGrid,1) || ...
            size(currentStack,2) ~= size(XGrid,2)
        error('Experiment %d: stack and grid dimensions do not match.', experimentNumber);
    end

    STD_TWO_ALL{experimentIndex} = double(currentStack);
    Grid_TWO_ALL{experimentIndex} = {XGrid, YGrid, xVector, yVector};
    CENTER1(experimentIndex,:) = Saved_CENTER1(savedIndex,:);
    CENTER2(experimentIndex,:) = Saved_CENTER2(savedIndex,:);
    LRadius1(experimentIndex) = Saved_LRadius1(savedIndex);
    LRadius2(experimentIndex) = Saved_LRadius2(savedIndex);
end

fprintf('\nLoaded two-spot experiments: %s\n', num2str(Experiments_Two'));

%% ============================================================
% 3. NORMALIZE TWO-SPOT CUMULATIVE STD STACKS
%% ============================================================

Normalized_STD_TWO_ALL = cell(numberTwoSpotExperiments, 1);
numberFramesTwoSpot = zeros(numberTwoSpotExperiments, 1);

for experimentIndex = 1:numberTwoSpotExperiments
    currentStack = STD_TWO_ALL{experimentIndex};
    currentGrid = Grid_TWO_ALL{experimentIndex};
    XGrid = currentGrid{1};
    YGrid = currentGrid{2};

    distances1 = hypot(XGrid-CENTER1(experimentIndex,1), ...
                       YGrid-CENTER1(experimentIndex,2));
    distances2 = hypot(XGrid-CENTER2(experimentIndex,1), ...
                       YGrid-CENTER2(experimentIndex,2));

    normalizationMask = ...
        distances1 <= LRadius1(experimentIndex) | ...
        distances2 <= LRadius2(experimentIndex);

    Normalized_STD_TWO_ALL{experimentIndex} = ...
        percentileNormalizeStack(currentStack, normalizationMask, ...
        normalizationLowerPercentile, normalizationUpperPercentile);

    numberFramesTwoSpot(experimentIndex) = size(currentStack, 3);
end

%% ============================================================
% 4. PROCESS SINGLE-SPOT TIFFS USING ORIGINAL CUMULATIVE METHOD
%% ============================================================

fprintf('\n============================================\n');
fprintf('PROCESSING SINGLE-SPOT REFERENCE EXPERIMENTS\n');
fprintf('============================================\n');

Experiments_Single = double(Experiments_Single(:));
numberSingleSpotExperiments = numel(Experiments_Single);

Normalized_STD_ONE_ALL = cell(numberSingleSpotExperiments, 1);
Grid_ONE_ALL = cell(numberSingleSpotExperiments, 1);
Center_ONE_ALL = nan(numberSingleSpotExperiments, 2);
LRadius_One = nan(numberSingleSpotExperiments, 1);
numberFramesSingleSpot = zeros(numberSingleSpotExperiments, 1);
SingleConditionPassed = false(numberSingleSpotExperiments, 1);
SingleConditionMessage = strings(numberSingleSpotExperiments, 1);
SingleTIFFFile = strings(numberSingleSpotExperiments, 1);
SingleMetadataFile = strings(numberSingleSpotExperiments, 1);

for experimentIndex = 1:numberSingleSpotExperiments
    experimentNumber = Experiments_Single(experimentIndex);
    fprintf('\nSingle-spot experiment %d\n', experimentNumber);

    experimentDirectory = fullfile(singleDataBasePath, ...
        [singleFolderPrefix num2str(experimentNumber)]);

    if ~isfolder(experimentDirectory)
        SingleConditionMessage(experimentIndex) = "Experiment folder missing";
        warning('Folder missing: %s', experimentDirectory);
        continue;
    end

    metadataFiles = removeInvalidFiles(dir(fullfile(experimentDirectory, '*.txt')));
    metadataParsed = false;

    for metadataIndex = 1:numel(metadataFiles)
        metadataPath = fullfile(experimentDirectory, metadataFiles(metadataIndex).name);
        try
            [center, xUnit_um, yUnit_um, laserRadius_um] = ...
                parseMetadata(metadataPath, 'single');
            metadataParsed = true;
            SingleMetadataFile(experimentIndex) = string(metadataFiles(metadataIndex).name);
            break;
        catch ME
            fprintf('  Rejected metadata "%s": %s\n', ...
                metadataFiles(metadataIndex).name, ME.message);
        end
    end

    if ~metadataParsed
        SingleConditionMessage(experimentIndex) = "Metadata parsing failed";
        continue;
    end

    tiffFiles = findTIFFCandidates(experimentDirectory);
    validTIFFFound = false;
    tiffPath = '';
    tiffInfo = [];

    for candidateIndex = 1:numel(tiffFiles)
        candidatePath = fullfile(experimentDirectory, tiffFiles(candidateIndex).name);
        try
            candidateInfo = imfinfo(candidatePath);
            testImage = imread(candidatePath, 1);
            if isempty(candidateInfo) || isempty(testImage) || ndims(testImage) ~= 2
                continue;
            end
            tiffPath = candidatePath;
            tiffInfo = candidateInfo;
            validTIFFFound = true;
            SingleTIFFFile(experimentIndex) = string(tiffFiles(candidateIndex).name);
            fprintf('  Selected TIFF: %s\n', tiffFiles(candidateIndex).name);
            break;
        catch ME
            fprintf('  Rejected TIFF "%s": %s\n', ...
                tiffFiles(candidateIndex).name, ME.message);
        end
    end

    if ~validTIFFFound
        SingleConditionMessage(experimentIndex) = "No readable TIFF";
        continue;
    end

    availableFrames = numel(tiffInfo);
    firstImage = double(imread(tiffPath, 1));
    [imageHeight, imageWidth] = size(firstImage);

    xVector = (1:imageWidth)*xUnit_um;
    yVector = (1:imageHeight)*yUnit_um;
    [XGrid, YGrid] = meshgrid(xVector, yVector);

    rawIntensityStack = zeros(imageHeight, imageWidth, availableFrames, 'double');
    cumulativeSTDStack = zeros(imageHeight, imageWidth, availableFrames, 'double');

    for frameIndex = 1:availableFrames
        rawIntensityStack(:,:,frameIndex) = double(imread(tiffPath, frameIndex));
        cumulativeSTDStack(:,:,frameIndex) = ...
            std(rawIntensityStack(:,:,1:frameIndex), 0, 3);
    end

    referenceNormalizationMask = ...
        hypot(XGrid-center(1), YGrid-center(2)) <= laserRadius_um;

    normalizedSTDStack = percentileNormalizeStack( ...
        cumulativeSTDStack, referenceNormalizationMask, ...
        normalizationLowerPercentile, normalizationUpperPercentile);

    Normalized_STD_ONE_ALL{experimentIndex} = normalizedSTDStack;
    Grid_ONE_ALL{experimentIndex} = {XGrid, YGrid, xVector, yVector};
    Center_ONE_ALL(experimentIndex,:) = center;
    LRadius_One(experimentIndex) = laserRadius_um;
    numberFramesSingleSpot(experimentIndex) = availableFrames;
    SingleConditionPassed(experimentIndex) = true;
    SingleConditionMessage(experimentIndex) = "Passed";

    fprintf('  Loaded %d cumulative frames; R = %.3f um\n', ...
        availableFrames, laserRadius_um);

    clear rawIntensityStack cumulativeSTDStack normalizedSTDStack;
end

if nnz(SingleConditionPassed) < 1
    error('No valid single-spot reference experiments were loaded.');
end

%% ============================================================
% 5. AVERAGE SINGLE-SPOT MAPS FRAME BY FRAME
%% ============================================================

validSingleIndices = find(SingleConditionPassed);
referenceSingleIndex = validSingleIndices(1);
referenceGrid = Grid_ONE_ALL{referenceSingleIndex};
referenceXVector = referenceGrid{3};
referenceYVector = referenceGrid{4};
referenceImageSize = size(Normalized_STD_ONE_ALL{referenceSingleIndex}(:,:,1));

for k = 1:numel(validSingleIndices)
    idx = validSingleIndices(k);
    if ~isequal(size(Normalized_STD_ONE_ALL{idx}(:,:,1)), referenceImageSize)
        error('Single-spot image sizes differ; old averaging method cannot be used.');
    end
end

maximumSingleAverageFrame = max(numberFramesSingleSpot(validSingleIndices));
Intensity_ONE_Frame = cell(1, maximumSingleAverageFrame);
NumberSingleExperimentsInAverage = zeros(1, maximumSingleAverageFrame);

for frameIndex = 1:maximumSingleAverageFrame
    accumulatedMap = zeros(referenceImageSize);
    validCount = 0;

    for k = 1:numel(validSingleIndices)
        idx = validSingleIndices(k);
        currentStack = Normalized_STD_ONE_ALL{idx};
        if frameIndex <= size(currentStack,3)
            accumulatedMap = accumulatedMap + currentStack(:,:,frameIndex);
            validCount = validCount + 1;
        end
    end

    if validCount > 0
        Intensity_ONE_Frame{frameIndex} = accumulatedMap/validCount;
    end
    NumberSingleExperimentsInAverage(frameIndex) = validCount;
end

avg_center_one = mean(Center_ONE_ALL(SingleConditionPassed,:), 1, 'omitnan');
avg_radius_one = mean(LRadius_One(SingleConditionPassed), 'omitnan');

%% ============================================================
% 6. DETERMINE FRAME RANGE
%% ============================================================

maximumAvailableFrame = min([numberFramesRequested, ...
    max(numberFramesTwoSpot), numel(Intensity_ONE_Frame)]);

if isempty(framesToAnalyze)
    framesToAnalyze = 1:maximumAvailableFrame;
else
    framesToAnalyze = unique(round(framesToAnalyze(:)'));
    framesToAnalyze = framesToAnalyze(framesToAnalyze >= 1 & ...
        framesToAnalyze <= maximumAvailableFrame);
end

if isempty(framesToAnalyze)
    error('No valid frames remain.');
end

numberAnalyzedFrames = numel(framesToAnalyze);
xTime_s = (framesToAnalyze-1)*frameInterval_s;

%% ============================================================
% 7. SELECT MANUAL CONTROL SQUARES
%% ============================================================

ControlRect1_um = nan(numberTwoSpotExperiments, 4);
ControlRect2_um = nan(numberTwoSpotExperiments, 4);
ControlSelectionCompleted = false(numberTwoSpotExperiments, 1);

if calculateControl && reuseSavedControlROIs && isfile(controlROIFile)
    savedControlData = load(controlROIFile);
    if isfield(savedControlData,'ControlRect1_um') && ...
            isfield(savedControlData,'ControlRect2_um') && ...
            isfield(savedControlData,'ControlSelectionCompleted') && ...
            isfield(savedControlData,'Experiments_Two') && ...
            isequal(double(savedControlData.Experiments_Two(:)), Experiments_Two(:))
        ControlRect1_um = savedControlData.ControlRect1_um;
        ControlRect2_um = savedControlData.ControlRect2_um;
        ControlSelectionCompleted = logical(savedControlData.ControlSelectionCompleted);
    end
end

if calculateControl
    for experimentIndex = 1:numberTwoSpotExperiments
        if ControlSelectionCompleted(experimentIndex)
            continue;
        end

        experimentNumber = Experiments_Two(experimentIndex);
        currentStack = Normalized_STD_TWO_ALL{experimentIndex};
        currentGrid = Grid_TWO_ALL{experimentIndex};
        XGrid = currentGrid{1};
        YGrid = currentGrid{2};

        availableDisplayFrames = framesToAnalyze(framesToAnalyze <= size(currentStack,3));
        if isempty(availableDisplayFrames)
            continue;
        end

        selectionFinished = false;
        while ~selectionFinished
            switch lower(controlSelectionDisplay)
                case 'temporalmean'
                    displayMap = mean(currentStack(:,:,availableDisplayFrames), 3, 'omitnan');
                    displayDescription = 'temporal-mean cumulative STD';
                case 'selectedframe'
                    selectedFrame = min(max(round(controlSelectionFrame),1), size(currentStack,3));
                    displayMap = currentStack(:,:,selectedFrame);
                    displayDescription = sprintf('frame %d', selectedFrame);
                otherwise
                    error('Invalid controlSelectionDisplay.');
            end

            figSelection = figure('Color','w','Position',[100 100 980 760]);
            axSelection = axes(figSelection);
            imagesc(axSelection, currentGrid{3}, currentGrid{4}, displayMap);
            axis(axSelection,'image');
            set(axSelection,'YDir','normal');
            colormap(axSelection, purpleProbabilityMap(256, rgb_Purple));
            colorbar(axSelection);
            hold(axSelection,'on');

            center1 = CENTER1(experimentIndex,:);
            center2 = CENTER2(experimentIndex,:);
            exclusionRadius1_um = controlExclusionMultiplier*LRadius1(experimentIndex);
            exclusionRadius2_um = controlExclusionMultiplier*LRadius2(experimentIndex);
            theta = linspace(0,2*pi,500);

            plot(axSelection, center1(1), center1(2), 'w+', 'MarkerSize',16,'LineWidth',2);
            plot(axSelection, center2(1), center2(2), 'w+', 'MarkerSize',16,'LineWidth',2);
            plot(axSelection, center1(1)+exclusionRadius1_um*cos(theta), ...
                center1(2)+exclusionRadius1_um*sin(theta), 'w--','LineWidth',1.5);
            plot(axSelection, center2(1)+exclusionRadius2_um*cos(theta), ...
                center2(2)+exclusionRadius2_um*sin(theta), 'w--','LineWidth',1.5);

            title(axSelection, {sprintf('Experiment %d: %s', ...
                experimentNumber, displayDescription), ...
                'Draw two equal nonactivated squares'}, 'FontWeight','normal');
            xlabel(axSelection,'X (\mum)');
            ylabel(axSelection,'Y (\mum)');
            formatPaperAxes(axSelection,FZ,rgb_Black);

            defaultPosition1 = chooseDefaultSquarePosition(XGrid,YGrid,fixedControlSquareSide_um);
            square1 = drawrectangle(axSelection,'Position',defaultPosition1, ...
                'Color',rgb_Cyan,'LineWidth',2,'FixedAspectRatio',true);
            wait(square1);
            position1 = square1.Position;

            if isempty(fixedControlSquareSide_um)
                sideLength_um = min(position1(3),position1(4));
            else
                sideLength_um = fixedControlSquareSide_um;
            end

            position1(3:4) = sideLength_um;
            position1 = constrainRectangleToImage(position1,XGrid,YGrid);
            square1.Position = position1;

            defaultPosition2 = position1;
            defaultPosition2(1) = position1(1)+1.5*sideLength_um;
            defaultPosition2 = constrainRectangleToImage(defaultPosition2,XGrid,YGrid);

            square2 = drawrectangle(axSelection,'Position',defaultPosition2, ...
                'Color',rgb_Orange,'LineWidth',2,'FixedAspectRatio',true);
            wait(square2);
            position2 = square2.Position;
            position2(3:4) = sideLength_um;
            position2 = constrainRectangleToImage(position2,XGrid,YGrid);
            square2.Position = position2;

            warningText = '';
            if rectangleIntersectsCircle(position1,center1,exclusionRadius1_um) || ...
                    rectangleIntersectsCircle(position1,center2,exclusionRadius2_um) || ...
                    rectangleIntersectsCircle(position2,center1,exclusionRadius1_um) || ...
                    rectangleIntersectsCircle(position2,center2,exclusionRadius2_um)
                warningText = [warningText sprintf('\nA square overlaps an exclusion circle.')];
            end
            if rectanglesOverlap(position1,position2)
                warningText = [warningText sprintf('\nThe two squares overlap.')];
            end
            if isempty(warningText)
                warningText = '\nNo geometric warning was detected.';
            end

            choice = questdlg(sprintf(['Experiment %d\nSquare side: %.2f um\n%s\n\n', ...
                'Use these regions?'], experimentNumber, sideLength_um, warningText), ...
                'Confirm control regions','Use regions','Redraw','Skip experiment','Use regions');

            switch choice
                case 'Use regions'
                    ControlRect1_um(experimentIndex,:) = position1;
                    ControlRect2_um(experimentIndex,:) = position2;
                    ControlSelectionCompleted(experimentIndex) = true;
                    selectionFinished = true;
                case 'Redraw'
                    selectionFinished = false;
                otherwise
                    selectionFinished = true;
            end

            if isgraphics(figSelection)
                close(figSelection);
            end
        end
    end

    save(controlROIFile,'ControlRect1_um','ControlRect2_um', ...
        'ControlSelectionCompleted','Experiments_Two');
end

%% ============================================================
% 8. PREALLOCATE KL RESULTS
%% ============================================================

Activated_KL_Matrix = nan(numberTwoSpotExperiments, maximumAvailableFrame);
Control_KL_Matrix = nan(numberTwoSpotExperiments, maximumAvailableFrame);
Activated_KL_Q1 = nan(size(Activated_KL_Matrix));
Activated_KL_Q2 = nan(size(Activated_KL_Matrix));
Control_KL_C1 = nan(size(Control_KL_Matrix));
Control_KL_C2 = nan(size(Control_KL_Matrix));
Activated_NumberSamples_Q1 = nan(size(Activated_KL_Matrix));
Activated_NumberSamples_Q2 = nan(size(Activated_KL_Matrix));
Control_NumberSamples_C1 = nan(size(Control_KL_Matrix));
Control_NumberSamples_C2 = nan(size(Control_KL_Matrix));
Activated_KL_Frame = cell(numberTwoSpotExperiments, maximumAvailableFrame);
Control_KL_Frame = cell(numberTwoSpotExperiments, maximumAvailableFrame);
ReferenceSamplesByExperimentFrame = cell(numberTwoSpotExperiments, maximumAvailableFrame);
ReferencePDF_ALL = cell(numberTwoSpotExperiments, maximumAvailableFrame);
ActivatedPDF_Q1_ALL = cell(size(Activated_KL_Frame));
ActivatedPDF_Q2_ALL = cell(size(Activated_KL_Frame));
ControlPDF_C1_ALL = cell(size(Control_KL_Frame));
ControlPDF_C2_ALL = cell(size(Control_KL_Frame));

%% ============================================================
% 9. CALCULATE ACTIVATED AND CONTROL KL
%% ============================================================

fprintf('\n============================================\n');
fprintf('CALCULATING ACTIVATED AND CONTROL KL\n');
fprintf('============================================\n');

for experimentIndex = 1:numberTwoSpotExperiments
    experimentNumber = Experiments_Two(experimentIndex);
    currentStack = Normalized_STD_TWO_ALL{experimentIndex};
    currentGrid = Grid_TWO_ALL{experimentIndex};
    X_two = currentGrid{1};
    Y_two = currentGrid{2};

    center1 = CENTER1(experimentIndex,:);
    center2 = CENTER2(experimentIndex,:);

    activatedMask1 = hypot(X_two-center1(1),Y_two-center1(2)) <= ...
        activatedROIRadiusMultiplier*LRadius1(experimentIndex);
    activatedMask2 = hypot(X_two-center2(1),Y_two-center2(2)) <= ...
        activatedROIRadiusMultiplier*LRadius2(experimentIndex);

    referenceMaskOnTwoGrid = ...
        hypot(X_two-avg_center_one(1),Y_two-avg_center_one(2)) <= ...
        singleReferenceRadiusMultiplier*avg_radius_one;

    if calculateControl && ControlSelectionCompleted(experimentIndex)
        controlMask1 = rectangleMaskFromPosition(X_two,Y_two,ControlRect1_um(experimentIndex,:));
        controlMask2 = rectangleMaskFromPosition(X_two,Y_two,ControlRect2_um(experimentIndex,:));
    else
        controlMask1 = false(size(X_two));
        controlMask2 = false(size(X_two));
    end

    validFrames = framesToAnalyze(framesToAnalyze <= size(currentStack,3));

    for frameIndex = validFrames
        if isempty(Intensity_ONE_Frame{frameIndex})
            continue;
        end

        F_interp = griddedInterpolant( ...
            {referenceYVector,referenceXVector}, ...
            Intensity_ONE_Frame{frameIndex}, ...
            'linear','none');

        intensity_matrix_one_interp = F_interp(Y_two,X_two);
        referenceSamples = intensity_matrix_one_interp(referenceMaskOnTwoGrid);
        referenceSamples = referenceSamples(isfinite(referenceSamples));
        ReferenceSamplesByExperimentFrame{experimentIndex,frameIndex} = referenceSamples;

        currentMap = currentStack(:,:,frameIndex);
        Q1Samples = currentMap(activatedMask1);
        Q2Samples = currentMap(activatedMask2);
        Q1Samples = Q1Samples(isfinite(Q1Samples));
        Q2Samples = Q2Samples(isfinite(Q2Samples));

        if removeOutliersBeforeKL
            referenceSamples = removeOutliersIQR(referenceSamples,outlierIQRMultiplier);
            Q1Samples = removeOutliersIQR(Q1Samples,outlierIQRMultiplier);
            Q2Samples = removeOutliersIQR(Q2Samples,outlierIQRMultiplier);
        end

        Activated_NumberSamples_Q1(experimentIndex,frameIndex) = numel(Q1Samples);
        Activated_NumberSamples_Q2(experimentIndex,frameIndex) = numel(Q2Samples);

        if frameIndex == 1
            Activated_KL_Q1(experimentIndex,frameIndex) = 0;
            Activated_KL_Q2(experimentIndex,frameIndex) = 0;
            Activated_KL_Matrix(experimentIndex,frameIndex) = 0;
            Activated_KL_Frame{experimentIndex,frameIndex} = 0;
        else
            [kl1,kl2,P,Q1,Q2] = compareTwoDistributionsOriginalMethod( ...
                referenceSamples,Q1Samples,Q2Samples,evaluationPoints, ...
                kdeEpsilon,minimumSamplesPerDistribution);

            Activated_KL_Q1(experimentIndex,frameIndex) = kl1;
            Activated_KL_Q2(experimentIndex,frameIndex) = kl2;

            if isfinite(kl1) && isfinite(kl2)
                Activated_KL_Matrix(experimentIndex,frameIndex) = mean([kl1 kl2]);
                Activated_KL_Frame{experimentIndex,frameIndex} = ...
                    Activated_KL_Matrix(experimentIndex,frameIndex);
            end

            ReferencePDF_ALL{experimentIndex,frameIndex} = P;
            ActivatedPDF_Q1_ALL{experimentIndex,frameIndex} = Q1;
            ActivatedPDF_Q2_ALL{experimentIndex,frameIndex} = Q2;
        end

        if calculateControl && ControlSelectionCompleted(experimentIndex)
            C1Samples = currentMap(controlMask1);
            C2Samples = currentMap(controlMask2);
            C1Samples = C1Samples(isfinite(C1Samples));
            C2Samples = C2Samples(isfinite(C2Samples));

            if removeOutliersBeforeKL
                C1Samples = removeOutliersIQR(C1Samples,outlierIQRMultiplier);
                C2Samples = removeOutliersIQR(C2Samples,outlierIQRMultiplier);
            end

            Control_NumberSamples_C1(experimentIndex,frameIndex) = numel(C1Samples);
            Control_NumberSamples_C2(experimentIndex,frameIndex) = numel(C2Samples);

            if frameIndex == 1
                Control_KL_C1(experimentIndex,frameIndex) = 0;
                Control_KL_C2(experimentIndex,frameIndex) = 0;
                Control_KL_Matrix(experimentIndex,frameIndex) = 0;
                Control_KL_Frame{experimentIndex,frameIndex} = 0;
            else
                [klc1,klc2,~,C1,C2] = compareTwoDistributionsOriginalMethod( ...
                    referenceSamples,C1Samples,C2Samples,evaluationPoints, ...
                    kdeEpsilon,minimumSamplesPerDistribution);

                Control_KL_C1(experimentIndex,frameIndex) = klc1;
                Control_KL_C2(experimentIndex,frameIndex) = klc2;

                if isfinite(klc1) && isfinite(klc2)
                    Control_KL_Matrix(experimentIndex,frameIndex) = mean([klc1 klc2]);
                    Control_KL_Frame{experimentIndex,frameIndex} = ...
                        Control_KL_Matrix(experimentIndex,frameIndex);
                end

                ControlPDF_C1_ALL{experimentIndex,frameIndex} = C1;
                ControlPDF_C2_ALL{experimentIndex,frameIndex} = C2;
            end
        end

        fprintf('Exp %d frame %d: activated KL = %.6f, control KL = %.6f\n', ...
            experimentNumber,frameIndex, ...
            Activated_KL_Matrix(experimentIndex,frameIndex), ...
            Control_KL_Matrix(experimentIndex,frameIndex));
    end
end

%% ============================================================
% 10. FRAMEWISE STATISTICS
%% ============================================================

[ActivatedMeanKL,ActivatedSDKL,ActivatedSEMKL,ActivatedN] = ...
    calculateFrameStatistics(Activated_KL_Matrix,framesToAnalyze);

[ControlMeanKL,ControlSDKL,ControlSEMKL,ControlN] = ...
    calculateFrameStatistics(Control_KL_Matrix,framesToAnalyze);

%% ============================================================
% 11. FIGURES
%% ============================================================

figActivatedKL = plotEvolutionFigure(xTime_s, ...
    Activated_KL_Matrix(:,framesToAnalyze),ActivatedMeanKL,ActivatedSDKL, ...
    rgb_Purple,rgb_LightPurple,rgb_ShadePurple,'-o', ...
    'Activated regions','KL divergence',FZ,rgb_Black,figurePosition);

figControlKL = plotEvolutionFigure(xTime_s, ...
    Control_KL_Matrix(:,framesToAnalyze),ControlMeanKL,ControlSDKL, ...
    rgb_Gray,rgb_LightGray,rgb_ShadeGray,'-s', ...
    'Nonactivated control regions','KL divergence',FZ,rgb_Black,figurePosition);

figKLComparison = figure('Color','w','Position',figurePosition);
ax = axes(figKLComparison);
hold(ax,'on');

activatedLower = max(ActivatedMeanKL-ActivatedSDKL,0);
activatedUpper = ActivatedMeanKL+ActivatedSDKL;
controlLower = max(ControlMeanKL-ControlSDKL,0);
controlUpper = ControlMeanKL+ControlSDKL;

valid = isfinite(activatedLower) & isfinite(activatedUpper);
if nnz(valid) >= 2
    fill(ax,[xTime_s(valid) fliplr(xTime_s(valid))], ...
        [activatedUpper(valid) fliplr(activatedLower(valid))], ...
        rgb_ShadePurple,'FaceAlpha',0.45,'EdgeColor','none', ...
        'HandleVisibility','off');
end

valid = isfinite(controlLower) & isfinite(controlUpper);
if nnz(valid) >= 2
    fill(ax,[xTime_s(valid) fliplr(xTime_s(valid))], ...
        [controlUpper(valid) fliplr(controlLower(valid))], ...
        rgb_ShadeGray,'FaceAlpha',0.45,'EdgeColor','none', ...
        'HandleVisibility','off');
end

plot(ax,xTime_s,ActivatedMeanKL,'-o','Color',rgb_Purple, ...
    'LineWidth',3.2,'MarkerSize',7,'MarkerFaceColor',rgb_Purple, ...
    'MarkerEdgeColor','none','DisplayName','Activated regions');

plot(ax,xTime_s,ControlMeanKL,'-s','Color',rgb_Gray, ...
    'LineWidth',2.8,'MarkerSize',7,'MarkerFaceColor',rgb_Gray, ...
    'MarkerEdgeColor','none','DisplayName','Nonactivated controls');

xlabel(ax,'Time (s)');
ylabel(ax,'KL divergence');
title(ax,'Activated and nonactivated regions','FontWeight','normal');
formatPaperAxes(ax,FZ,rgb_Black);
legend(ax,'Location','northwest','Box','off');
xlim(ax,[min(xTime_s) max(xTime_s)]);
hold(ax,'off');

%% ============================================================
% 12. BUILD TABLES
%% ============================================================

ExperimentColumn = [];
FrameColumn = [];
TimeColumn_s = [];
ActivatedKLColumn = [];
ControlKLColumn = [];

for experimentIndex = 1:numberTwoSpotExperiments
    for frameIndex = framesToAnalyze
        ExperimentColumn(end+1,1) = Experiments_Two(experimentIndex); %#ok<SAGROW>
        FrameColumn(end+1,1) = frameIndex; %#ok<SAGROW>
        TimeColumn_s(end+1,1) = (frameIndex-1)*frameInterval_s; %#ok<SAGROW>
        ActivatedKLColumn(end+1,1) = Activated_KL_Matrix(experimentIndex,frameIndex); %#ok<SAGROW>
        ControlKLColumn(end+1,1) = Control_KL_Matrix(experimentIndex,frameIndex); %#ok<SAGROW>
    end
end

FrameResolvedTable = table(ExperimentColumn,FrameColumn,TimeColumn_s, ...
    ActivatedKLColumn,ControlKLColumn, ...
    'VariableNames',{'Experiment','Frame','Time_s','ActivatedKL','ControlKL'});

FrameSummaryTable = table(framesToAnalyze(:),xTime_s(:), ...
    ActivatedMeanKL(:),ActivatedSDKL(:),ActivatedSEMKL(:),ActivatedN(:), ...
    ControlMeanKL(:),ControlSDKL(:),ControlSEMKL(:),ControlN(:), ...
    'VariableNames',{'Frame','Time_s','ActivatedMeanKL','ActivatedSDKL', ...
    'ActivatedSEMKL','ActivatedN','ControlMeanKL','ControlSDKL', ...
    'ControlSEMKL','ControlN'});

%% ============================================================
% 13. SAVE OUTPUTS
%% ============================================================

KL_Frame = Activated_KL_Frame;
controlKL_Frame = Control_KL_Frame;

save(fullfile(outputDir,'KL_Frame_ActivatedVsControl.mat'), ...
    'KL_Frame','controlKL_Frame','Activated_KL_Frame','Control_KL_Frame', ...
    'Activated_KL_Matrix','Control_KL_Matrix','Experiments_Two', ...
    'framesToAnalyze','frameInterval_s');

if saveExcel
    outputExcel = fullfile(outputDir,'Activated_vs_Control_KL_OriginalMethod.xlsx');
    if isfile(outputExcel), delete(outputExcel); end
    writetable(FrameResolvedTable,outputExcel,'Sheet','FrameResolved');
    writetable(FrameSummaryTable,outputExcel,'Sheet','FrameSummary');
    writematrix(Activated_KL_Matrix,outputExcel,'Sheet','ActivatedKL_Matrix');
    writematrix(Control_KL_Matrix,outputExcel,'Sheet','ControlKL_Matrix');
end

if saveMAT
    save(fullfile(outputDir,'Activated_vs_Control_KL_OriginalMethod.mat'), ...
        'Activated_KL_Matrix','Control_KL_Matrix', ...
        'Activated_KL_Frame','Control_KL_Frame','KL_Frame','controlKL_Frame', ...
        'Activated_KL_Q1','Activated_KL_Q2','Control_KL_C1','Control_KL_C2', ...
        'ControlRect1_um','ControlRect2_um','ControlSelectionCompleted', ...
        'Intensity_ONE_Frame','Normalized_STD_ONE_ALL','Normalized_STD_TWO_ALL', ...
        'Experiments_Single','Experiments_Two','framesToAnalyze','frameInterval_s', ...
        'Center_ONE_ALL','LRadius_One','avg_center_one','avg_radius_one', ...
        'CENTER1','CENTER2','LRadius1','LRadius2', ...
        'FrameResolvedTable','FrameSummaryTable','-v7.3');
end

if saveFigures
    saveVectorAndPNG(figActivatedKL,outputDir,'ActivatedRegions_KL_OriginalMethod');
    saveVectorAndPNG(figControlKL,outputDir,'NonactivatedControlRegions_KL_OriginalMethod');
    saveVectorAndPNG(figKLComparison,outputDir,'ActivatedVsControl_KL_OriginalMethod');
end

fprintf('\nAnalysis complete. Output folder:\n%s\n',outputDir);

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function [center,xUnit,yUnit,radius] = parseMetadata(filePath,spotType)
    if nargin < 2, spotType = 'single'; end
    fileID = fopen(filePath,'r');
    if fileID == -1, error('Cannot open metadata file: %s',filePath); end
    cleaner = onCleanup(@() fclose(fileID)); %#ok<NASGU>
    C = textscan(fileID,'%s','Delimiter','\n');
    lines = C{1};
    xUnit = parsePixelUnit(lines{5});
    yUnit = parsePixelUnit(lines{6});

    switch lower(spotType)
        case {'single','spot1'}, lineIndex = numel(lines)-4;
        case 'spot2', lineIndex = numel(lines)-5;
        otherwise, error('Invalid spotType.');
    end

    line = lines{lineIndex};
    s = strfind(line,'(');
    e = strfind(line,')');
    points = [];
    for k = 1:min(numel(s),numel(e))
        txt = line(s(k)+1:e(k)-1);
        p = strsplit(txt,',');
        if numel(p) == 2
            x = str2double(strtrim(p{1}));
            y = str2double(strtrim(p{2}));
            if isfinite(x) && isfinite(y)
                points(end+1,:) = [x y]; %#ok<AGROW>
            end
        end
    end

    if size(points,1) < 2, error('Could not parse activation coordinates.'); end
    center = [0.5*(points(1,1)+points(2,1))*xUnit, ...
              0.5*(points(1,2)+points(2,2))*yUnit];
    radius = (1/(2*sqrt(2)))*hypot((points(2,1)-points(1,1))*xUnit, ...
                                   (points(2,2)-points(1,2))*yUnit);
end

function unit = parsePixelUnit(line)
    p = strsplit(line,':');
    q = strsplit(strtrim(p{2}),'*');
    unit = str2double(strtrim(q{2}));
    if ~isfinite(unit) || unit <= 0, error('Invalid pixel unit.'); end
end

function normalizedStack = percentileNormalizeStack(rawStack,mask,pLow,pHigh)
    values = [];
    for t = 1:size(rawStack,3)
        m = rawStack(:,:,t);
        v = m(mask);
        values = [values; v(isfinite(v))]; %#ok<AGROW>
    end
    low = prctile(values,pLow);
    high = prctile(values,pHigh);
    range = high-low;
    if range <= eps
        normalizedStack = zeros(size(rawStack));
        return;
    end
    normalizedStack = (double(rawStack)-low)/range;
    normalizedStack(normalizedStack < 0) = 0;
    normalizedStack(normalizedStack > 1) = 1;
end

function [KL1,KL2,P,Q1,Q2] = compareTwoDistributionsOriginalMethod( ...
        referenceSamples,sampleSet1,sampleSet2,evaluationPoints,epsilonValue,minimumSamples)

    KL1 = NaN; KL2 = NaN; P = []; Q1 = []; Q2 = [];
    referenceSamples = referenceSamples(isfinite(referenceSamples));
    sampleSet1 = sampleSet1(isfinite(sampleSet1));
    sampleSet2 = sampleSet2(isfinite(sampleSet2));

    if numel(referenceSamples) < minimumSamples || ...
            numel(sampleSet1) < minimumSamples || ...
            numel(sampleSet2) < minimumSamples
        return;
    end

    Pz = zscoreSafeOriginal(referenceSamples);
    Q1z = zscoreSafeOriginal(sampleSet1);
    Q2z = zscoreSafeOriginal(sampleSet2);

    if isempty(Pz) || isempty(Q1z) || isempty(Q2z) || ...
            any(isnan(Pz)) || any(isnan(Q1z)) || any(isnan(Q2z))
        return;
    end

    P = ksdensity(Pz,evaluationPoints);  P = P/sum(P)+epsilonValue;
    Q1 = ksdensity(Q1z,evaluationPoints); Q1 = Q1/sum(Q1)+epsilonValue;
    Q2 = ksdensity(Q2z,evaluationPoints); Q2 = Q2/sum(Q2)+epsilonValue;

    KL1 = sum(Q1.*log(Q1./P),'all');
    KL2 = sum(Q2.*log(Q2./P),'all');
end

function z = zscoreSafeOriginal(x)
    s = std(x,'omitnan');
    if s > eps
        z = (x-mean(x,'omitnan'))/s;
    else
        z = NaN(size(x));
    end
end

function cleaned = removeOutliersIQR(x,multiplier)
    x = x(isfinite(x));
    q1 = prctile(x,25);
    q3 = prctile(x,75);
    iqrVal = q3-q1;
    cleaned = x(x >= q1-multiplier*iqrVal & x <= q3+multiplier*iqrVal);
end

function files = findTIFFCandidates(directoryPath)
    files = [dir(fullfile(directoryPath,'*.tif')); ...
             dir(fullfile(directoryPath,'*.tiff')); ...
             dir(fullfile(directoryPath,'*.TIF')); ...
             dir(fullfile(directoryPath,'*.TIFF'))];
    files = removeInvalidFiles(files);
    if isempty(files), return; end
    [~,order] = sort([files.bytes],'descend');
    files = files(order);
end

function files = removeInvalidFiles(files)
    if isempty(files), return; end
    files = files(~[files.isdir]);
    names = string({files.name});
    files = files(~startsWith(names,"._") & ~startsWith(names,"."));
    if isempty(files), return; end
    files = files([files.bytes] > 0);
end

function mask = rectangleMaskFromPosition(XGrid,YGrid,p)
    mask = XGrid >= p(1) & XGrid <= p(1)+p(3) & ...
           YGrid >= p(2) & YGrid <= p(2)+p(4);
end

function [XGrid,YGrid,xVector,yVector] = extractSpatialGrid(gridEntry)
    xInput = double(gridEntry{1});
    yInput = double(gridEntry{2});
    if isvector(xInput) && isvector(yInput)
        xVector = xInput(:)';
        yVector = yInput(:);
        [XGrid,YGrid] = meshgrid(xVector,yVector);
    else
        XGrid = xInput;
        YGrid = yInput;
        xVector = XGrid(1,:);
        yVector = YGrid(:,1);
    end
end

function p = chooseDefaultSquarePosition(XGrid,YGrid,side)
    xmin = min(XGrid,[],'all'); xmax = max(XGrid,[],'all');
    ymin = min(YGrid,[],'all'); ymax = max(YGrid,[],'all');
    if isempty(side), side = 0.15*min(xmax-xmin,ymax-ymin); end
    p = [xmin+0.08*(xmax-xmin), ymin+0.08*(ymax-ymin), side, side];
end

function p = constrainRectangleToImage(p,XGrid,YGrid)
    xmin = min(XGrid,[],'all'); xmax = max(XGrid,[],'all');
    ymin = min(YGrid,[],'all'); ymax = max(YGrid,[],'all');
    p(3) = min(p(3),xmax-xmin); p(4) = min(p(4),ymax-ymin);
    p(1) = min(max(p(1),xmin),xmax-p(3));
    p(2) = min(max(p(2),ymin),ymax-p(4));
end

function tf = rectangleIntersectsCircle(p,c,r)
    cx = min(max(c(1),p(1)),p(1)+p(3));
    cy = min(max(c(2),p(2)),p(2)+p(4));
    tf = (cx-c(1))^2 + (cy-c(2))^2 <= r^2;
end

function tf = rectanglesOverlap(a,b)
    tf = a(1) < b(1)+b(3) && a(1)+a(3) > b(1) && ...
         a(2) < b(2)+b(4) && a(2)+a(4) > b(2);
end

function [m,s,se,n] = calculateFrameStatistics(data,frames)
    m = nan(1,numel(frames)); s = m; se = m; n = zeros(1,numel(frames));
    for k = 1:numel(frames)
        v = data(:,frames(k));
        v = v(isfinite(v));
        n(k) = numel(v);
        if ~isempty(v), m(k) = mean(v); end
        if numel(v) >= 2
            s(k) = std(v,0);
            se(k) = s(k)/sqrt(numel(v));
        end
    end
end

function fig = plotEvolutionFigure(x,data,m,s,mainColor,individualColor,shadeColor, ...
        lineSpec,titleText,yLabelText,fontSize,blackColor,figurePosition)
    fig = figure('Color','w','Position',figurePosition);
    ax = axes(fig); hold(ax,'on');
    lower = max(m-s,0); upper = m+s;
    valid = isfinite(lower) & isfinite(upper);
    if nnz(valid) >= 2
        fill(ax,[x(valid) fliplr(x(valid))],[upper(valid) fliplr(lower(valid))], ...
            shadeColor,'FaceAlpha',0.55,'EdgeColor','none','DisplayName','Mean \pm SD');
    end
    for i = 1:size(data,1)
        plot(ax,x,data(i,:),lineSpec,'Color',individualColor,'LineWidth',1.4, ...
            'MarkerSize',4.5,'MarkerFaceColor',individualColor, ...
            'MarkerEdgeColor','none','HandleVisibility','off');
    end
    plot(ax,x,m,lineSpec,'Color',mainColor,'LineWidth',3.5,'MarkerSize',6.5, ...
        'MarkerFaceColor',mainColor,'MarkerEdgeColor','none','DisplayName','Mean');
    xlabel(ax,'Time (s)'); ylabel(ax,yLabelText);
    title(ax,titleText,'FontWeight','normal');
    formatPaperAxes(ax,fontSize,blackColor);
    legend(ax,'Location','northwest','Box','off');
    xlim(ax,[min(x) max(x)]);
    hold(ax,'off');
end

function formatPaperAxes(ax,fontSize,blackColor)
    set(ax,'FontName','Arial','FontSize',fontSize,'LineWidth',1.2, ...
        'TickDir','in','TickLength',[0.018 0.018],'Box','on', ...
        'XColor',blackColor,'YColor',blackColor,'Layer','top');
    grid(ax,'off');
end

function cmap = purpleProbabilityMap(n,mainPurple)
    anchors = [1 1 1; 0.75*[1 1 1]+0.25*mainPurple; mainPurple; 0.4*mainPurple];
    cmap = interp1(linspace(0,1,size(anchors,1)),anchors,linspace(0,1,n),'pchip');
    cmap = min(max(cmap,0),1);
end

function saveVectorAndPNG(fig,outputDir,baseName)
    if isempty(fig) || ~isgraphics(fig), return; end
    set(fig,'Renderer','painters');
    exportgraphics(fig,fullfile(outputDir,[baseName '.pdf']),'ContentType','vector');
    exportgraphics(fig,fullfile(outputDir,[baseName '.eps']),'ContentType','vector');
    exportgraphics(fig,fullfile(outputDir,[baseName '.png']),'Resolution',600);
end
