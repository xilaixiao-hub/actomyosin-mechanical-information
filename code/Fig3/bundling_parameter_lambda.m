%% ============================================================
% MULTI-EXPERIMENT BUNDLING PARAMETER OVER ALL FRAMES
% USING ONE MANUALLY DEFINED, SHARED ROI
%
% Bundling parameter definition:
%
%   I_avg(x,t) = mean_y[I(x,y,t)]
%
%                 std_x(I_avg(x,t))
%      BP(t) = -------------------------
%                mean_x(I_avg(x,t))
%
% All experiments use the same ROI and activation-circle geometry.
%% ============================================================

clearvars;
close all;
clc;

%% ============================================================
% 1. USER SETTINGS
%% ============================================================

basePath = ...
    'PATH_TO_DATA/Experiments/TwoLightSpots_Data&PIV';

experimentNumbers = [6 8 9 10 11 13];
experimentFolderPrefix = 'set';
tiffFileName = 'actin.tif';

frameInterval_s = 10;
framesToAnalyze = [];   % empty = all available frames

pixelSize_um = 0.160772;

% Shared manual ROI: [x, y, width, height] in pixels.
sharedROIPosition_px = [97.26706-20,92.44390+20,2*20.0965,172.82990-92.44390-2*20.0965]/pixelSize_um;

showActivationCircles = true;

% Shared activation geometry, [x,y] in micrometers.
activationCenters_um = [ ...
    97.26706, 92.44390; ...
    97.26706,172.82990];

activationRadii_um = [20.0965;20.0965];

% Optional experiment-specific shifts: [deltaX_px,deltaY_px].
% Leave all zeros when the image locations are identical.
experimentShift_px = zeros(numel(experimentNumbers),2);

% Diagnostic image.
diagnosticExperimentIndex = 1;
diagnosticFrame = 1;
displayPercentiles = [1,99];

% Group plot.
showIndividualTraces = true;
errorRegion = 'SD';   % 'SD' or 'SEM'

FZ = 18;

rgbMain = [105,55,130]/255;
rgbIndividual = [202,178,214]/255;
rgbShade = [225,211,232]/255;
rgbRectangle = [255,215,70]/255;
rgbCircle = [1,1,1];
rgbBlack = [30,30,30]/255;

figurePosition = [100,100,900,680];

outputDir = fullfile(pwd,'MultiExperiment_BundlingParameter_AllFrames');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

saveMAT = true;
saveExcel = true;
saveFigures = true;

%% ============================================================
% 2. VALIDATE SETTINGS
%% ============================================================

experimentNumbers = double(experimentNumbers(:));
numberExperiments = numel(experimentNumbers);

if size(experimentShift_px,1) ~= numberExperiments || ...
        size(experimentShift_px,2) ~= 2
    error('experimentShift_px must contain one [dx,dy] row per experiment.');
end

if numel(sharedROIPosition_px) ~= 4
    error('sharedROIPosition_px must be [x,y,width,height].');
end

sharedROIPosition_px = double(sharedROIPosition_px(:)');

if any(~isfinite(sharedROIPosition_px)) || ...
        sharedROIPosition_px(3) <= 0 || sharedROIPosition_px(4) <= 0
    error('The shared ROI is invalid.');
end

activationCenters_px = activationCenters_um/pixelSize_um;
activationRadii_px = activationRadii_um/pixelSize_um;

%% ============================================================
% 3. INSPECT EXPERIMENT FILES
%% ============================================================

TIFFPaths = strings(numberExperiments,1);
NumberAvailableFrames = zeros(numberExperiments,1);
ConditionPassed = false(numberExperiments,1);
ConditionMessage = strings(numberExperiments,1);

ImageRows = nan(numberExperiments,1);
ImageColumns = nan(numberExperiments,1);

fprintf('\n============================================\n');
fprintf('CHECKING EXPERIMENT FILES\n');
fprintf('============================================\n');

for e = 1:numberExperiments

    experimentNumber = experimentNumbers(e);

    experimentDirectory = fullfile( ...
        basePath,[experimentFolderPrefix num2str(experimentNumber)]);

    tiffPath = fullfile(experimentDirectory,tiffFileName);
    TIFFPaths(e) = string(tiffPath);

    fprintf('\nExperiment %d\n',experimentNumber);

    if ~isfile(tiffPath)
        ConditionMessage(e) = "TIFF file missing";
        warning('TIFF file not found:\n%s',tiffPath);
        continue;
    end

    try
        imageInfo = imfinfo(tiffPath);
        firstImage = double(imread(tiffPath,1,'Info',imageInfo));

        if ndims(firstImage) == 3
            firstImage = mean(firstImage,3);
        end

        [numberRows,numberColumns] = size(firstImage);

        currentROIPosition = sharedROIPosition_px + [ ...
            experimentShift_px(e,1),experimentShift_px(e,2),0,0];

        [xStart,yStart,xEnd,yEnd] = rectangleToPixelBounds( ...
            currentROIPosition,numberColumns,numberRows);

        ImageRows(e) = numberRows;
        ImageColumns(e) = numberColumns;
        NumberAvailableFrames(e) = numel(imageInfo);

        ConditionPassed(e) = true;
        ConditionMessage(e) = "Passed";

        fprintf('  image size: %d rows x %d columns\n',numberRows,numberColumns);
        fprintf('  available frames: %d\n',NumberAvailableFrames(e));
        fprintf('  ROI bounds: x = %d:%d, y = %d:%d\n', ...
            xStart,xEnd,yStart,yEnd);

    catch ME
        ConditionMessage(e) = string(ME.message);
        warning('Experiment %d failed: %s',experimentNumber,ME.message);
    end
end

if nnz(ConditionPassed) < 1
    error('No valid experiments were found.');
end

%% ============================================================
% 4. DETERMINE FRAME RANGE
%% ============================================================

maximumFrame = max(NumberAvailableFrames);

if isempty(framesToAnalyze)
    framesToAnalyze = 1:maximumFrame;
else
    framesToAnalyze = unique(round(double(framesToAnalyze(:)')));
    framesToAnalyze = framesToAnalyze( ...
        framesToAnalyze >= 1 & framesToAnalyze <= maximumFrame);
end

if isempty(framesToAnalyze)
    error('No valid frames remain.');
end

time_s = (framesToAnalyze-1)*frameInterval_s;

%% ============================================================
% 5. PREALLOCATE RESULTS
%% ============================================================

BP_Matrix = nan(numberExperiments,maximumFrame);
MeanIntensity_Matrix = nan(numberExperiments,maximumFrame);
STDIntensity_Matrix = nan(numberExperiments,maximumFrame);

NumberHorizontalScans = nan(numberExperiments,1);
NumberPointsPerScan = nan(numberExperiments,1);

ROIPositionUsed_px = nan(numberExperiments,4);
ROIBounds_px = nan(numberExperiments,4);

AverageProfile_ALL = cell(numberExperiments,maximumFrame);
STDProfileAcrossRows_ALL = cell(numberExperiments,maximumFrame);

%% ============================================================
% 6. CALCULATE BP FOR EVERY FRAME
%% ============================================================

fprintf('\n============================================\n');
fprintf('CALCULATING BUNDLING PARAMETER\n');
fprintf('============================================\n');

for e = 1:numberExperiments

    if ~ConditionPassed(e)
        continue;
    end

    experimentNumber = experimentNumbers(e);
    tiffPath = char(TIFFPaths(e));
    imageInfo = imfinfo(tiffPath);
    numberFramesCurrent = numel(imageInfo);

    fprintf('\nExperiment %d\n',experimentNumber);

    currentROIPosition = sharedROIPosition_px + [ ...
        experimentShift_px(e,1),experimentShift_px(e,2),0,0];

    ROIPositionUsed_px(e,:) = currentROIPosition;

    firstImage = double(imread(tiffPath,1,'Info',imageInfo));
    if ndims(firstImage) == 3
        firstImage = mean(firstImage,3);
    end

    [numberRows,numberColumns] = size(firstImage);

    [xStart,yStart,xEnd,yEnd] = rectangleToPixelBounds( ...
        currentROIPosition,numberColumns,numberRows);

    ROIBounds_px(e,:) = [xStart,xEnd,yStart,yEnd];
    NumberHorizontalScans(e) = yEnd-yStart+1;
    NumberPointsPerScan(e) = xEnd-xStart+1;

    validFrames = framesToAnalyze(framesToAnalyze <= numberFramesCurrent);

    for frameIndex = validFrames

        currentImage = double(imread( ...
            tiffPath,frameIndex,'Info',imageInfo));

        if ndims(currentImage) == 3
            currentImage = mean(currentImage,3);
        end

        currentROI = currentImage(yStart:yEnd,xStart:xEnd);

        averageProfile = mean(currentROI,1,'omitnan');
        stdProfileAcrossRows = std(currentROI,0,1,'omitnan');

        profileMean = mean(averageProfile,'omitnan');
        profileSTD = std(averageProfile,0,'omitnan');

        if isfinite(profileMean) && profileMean > 0
            BP_Matrix(e,frameIndex) = profileSTD/profileMean;
        end

        MeanIntensity_Matrix(e,frameIndex) = profileMean;
        STDIntensity_Matrix(e,frameIndex) = profileSTD;

        AverageProfile_ALL{e,frameIndex} = averageProfile;
        STDProfileAcrossRows_ALL{e,frameIndex} = stdProfileAcrossRows;

        fprintf('  Frame %d: mean = %.6g, std = %.6g, BP = %.6f\n', ...
            frameIndex,profileMean,profileSTD,BP_Matrix(e,frameIndex));
    end
end

%% ============================================================
% 7. GROUP STATISTICS
%% ============================================================

MeanBP = nan(1,maximumFrame);
SDBP = nan(1,maximumFrame);
SEMBP = nan(1,maximumFrame);
NumberValidExperimentsByFrame = zeros(1,maximumFrame);

for frameIndex = framesToAnalyze

    values = BP_Matrix(:,frameIndex);
    values = values(isfinite(values));

    NumberValidExperimentsByFrame(frameIndex) = numel(values);

    if ~isempty(values)
        MeanBP(frameIndex) = mean(values);
    end

    if numel(values) >= 2
        SDBP(frameIndex) = std(values,0);
        SEMBP(frameIndex) = SDBP(frameIndex)/sqrt(numel(values));
    end
end

%% ============================================================
% 8. BP EVOLUTION FIGURE
%% ============================================================

plotData = BP_Matrix(:,framesToAnalyze);
meanPlot = MeanBP(framesToAnalyze);
sdPlot = SDBP(framesToAnalyze);
semPlot = SEMBP(framesToAnalyze);

switch upper(errorRegion)
    case 'SD'
        errorValues = sdPlot;
        errorLabel = 'Mean \pm SD';
    case 'SEM'
        errorValues = semPlot;
        errorLabel = 'Mean \pm SEM';
    otherwise
        error('errorRegion must be ''SD'' or ''SEM''.');
end

figBP = figure('Color','w','Position',figurePosition);
ax = axes(figBP);
hold(ax,'on');

validShade = isfinite(meanPlot) & isfinite(errorValues);

if nnz(validShade) >= 2
    lowerBound = max(meanPlot(validShade)-errorValues(validShade),0);
    upperBound = meanPlot(validShade)+errorValues(validShade);

    fill(ax, ...
        [time_s(validShade),fliplr(time_s(validShade))], ...
        [upperBound,fliplr(lowerBound)], ...
        rgbShade,'FaceAlpha',0.55,'EdgeColor','none', ...
        'DisplayName',errorLabel);
end

if showIndividualTraces
    for e = 1:numberExperiments
        plot(ax,time_s,plotData(e,:),'-o', ...
            'Color',rgbIndividual,'LineWidth',1.3, ...
            'MarkerSize',4,'MarkerFaceColor',rgbIndividual, ...
            'MarkerEdgeColor','none','HandleVisibility','off');
    end
end

plot(ax,time_s,meanPlot,'-o', ...
    'Color',rgbMain,'LineWidth',3.4, ...
    'MarkerSize',6.5,'MarkerFaceColor',rgbMain, ...
    'MarkerEdgeColor','none','DisplayName','Mean');

xlabel(ax,'Time (s)');
ylabel(ax,'Bundling parameter, std(I_{avg})/mean(I_{avg})');
title(ax,'Averaged-line-scan bundling parameter','FontWeight','normal');
legend(ax,'Location','best','Box','off');
formatAxes(ax,FZ,rgbBlack);

if numel(time_s) >= 2
    xlim(ax,[min(time_s),max(time_s)]);
end

hold(ax,'off');

%% ============================================================
% 9. DIAGNOSTIC FULL-FIELD FIGURE
%% ============================================================

figDiagnostic = [];

if diagnosticExperimentIndex >= 1 && ...
        diagnosticExperimentIndex <= numberExperiments && ...
        ConditionPassed(diagnosticExperimentIndex)

    e = diagnosticExperimentIndex;
    tiffPath = char(TIFFPaths(e));
    imageInfo = imfinfo(tiffPath);

    selectedDiagnosticFrame = min( ...
        max(round(diagnosticFrame),1),numel(imageInfo));

    diagnosticImage = double(imread( ...
        tiffPath,selectedDiagnosticFrame,'Info',imageInfo));

    if ndims(diagnosticImage) == 3
        diagnosticImage = mean(diagnosticImage,3);
    end

    shiftedActivationCenters_px = activationCenters_px + ...
        experimentShift_px(e,:);

    figDiagnostic = figure('Color','w','Position',[100,100,760,700]);
    axDiagnostic = axes(figDiagnostic);

    displayPurpleImage(axDiagnostic,diagnosticImage,displayPercentiles);
    hold(axDiagnostic,'on');

    if showActivationCircles
        drawActivationCircles(axDiagnostic, ...
            shiftedActivationCenters_px,activationRadii_px, ...
            rgbCircle,2.5);
    end

    rectangle(axDiagnostic, ...
        'Position',ROIPositionUsed_px(e,:), ...
        'EdgeColor',rgbRectangle,'LineWidth',3);

    hold(axDiagnostic,'off');

    xlabel(axDiagnostic,'x (pixels)');
    ylabel(axDiagnostic,'y (pixels)');
    title(axDiagnostic,sprintf( ...
        'Experiment %d, frame %d: shared ROI', ...
        experimentNumbers(e),selectedDiagnosticFrame), ...
        'FontWeight','normal');

    formatAxes(axDiagnostic,FZ,rgbBlack);
end

%% ============================================================
% 10. HEATMAP
%% ============================================================

figHeatmap = figure('Color','w','Position',[100,100,960,580]);
axHeatmap = axes(figHeatmap);

imagesc(axHeatmap,time_s,1:numberExperiments, ...
    BP_Matrix(:,framesToAnalyze));

set(axHeatmap,'YDir','normal');
colormap(axHeatmap,purpleColormap(256));

cb = colorbar(axHeatmap);
cb.Label.String = 'Bundling parameter';

xlabel(axHeatmap,'Time (s)');
ylabel(axHeatmap,'Experiment');
yticks(axHeatmap,1:numberExperiments);
yticklabels(axHeatmap,string(experimentNumbers));
title(axHeatmap,'Bundling parameter by experiment and time', ...
    'FontWeight','normal');
formatAxes(axHeatmap,FZ,rgbBlack);

%% ============================================================
% 11. BUILD OUTPUT TABLES
%% ============================================================

ExperimentColumn = [];
FrameColumn = [];
TimeColumn_s = [];
BPColumn = [];
MeanIntensityColumn = [];
STDIntensityColumn = [];

for e = 1:numberExperiments
    for frameIndex = framesToAnalyze
        ExperimentColumn(end+1,1) = experimentNumbers(e); %#ok<SAGROW>
        FrameColumn(end+1,1) = frameIndex; %#ok<SAGROW>
        TimeColumn_s(end+1,1) = (frameIndex-1)*frameInterval_s; %#ok<SAGROW>
        BPColumn(end+1,1) = BP_Matrix(e,frameIndex); %#ok<SAGROW>
        MeanIntensityColumn(end+1,1) = MeanIntensity_Matrix(e,frameIndex); %#ok<SAGROW>
        STDIntensityColumn(end+1,1) = STDIntensity_Matrix(e,frameIndex); %#ok<SAGROW>
    end
end

FrameResolvedTable = table( ...
    ExperimentColumn,FrameColumn,TimeColumn_s,BPColumn, ...
    MeanIntensityColumn,STDIntensityColumn, ...
    'VariableNames',{ ...
        'Experiment','Frame','Time_s','BundlingParameter', ...
        'MeanOfAveragedProfile','STDOfAveragedProfile'});

FrameSummaryTable = table( ...
    framesToAnalyze(:),time_s(:), ...
    MeanBP(framesToAnalyze)',SDBP(framesToAnalyze)', ...
    SEMBP(framesToAnalyze)', ...
    NumberValidExperimentsByFrame(framesToAnalyze)', ...
    'VariableNames',{ ...
        'Frame','Time_s','MeanBundlingParameter', ...
        'SDBundlingParameter','SEMBundlingParameter', ...
        'NumberValidExperiments'});

ExperimentSettingsTable = table( ...
    experimentNumbers,TIFFPaths,NumberAvailableFrames, ...
    ROIPositionUsed_px(:,1),ROIPositionUsed_px(:,2), ...
    ROIPositionUsed_px(:,3),ROIPositionUsed_px(:,4), ...
    ROIBounds_px(:,1),ROIBounds_px(:,2), ...
    ROIBounds_px(:,3),ROIBounds_px(:,4), ...
    NumberHorizontalScans,NumberPointsPerScan, ...
    ConditionPassed,ConditionMessage, ...
    'VariableNames',{ ...
        'Experiment','TIFFPath','NumberAvailableFrames', ...
        'ROI_X_px','ROI_Y_px','ROI_Width_px','ROI_Height_px', ...
        'XStart_px','XEnd_px','YStart_px','YEnd_px', ...
        'NumberHorizontalScans','NumberPointsPerScan', ...
        'ConditionPassed','ConditionMessage'});

%% ============================================================
% 12. SAVE OUTPUTS
%% ============================================================

if saveExcel

    outputExcel = fullfile(outputDir, ...
        'MultiExperiment_BundlingParameter_AllFrames.xlsx');

    if isfile(outputExcel)
        delete(outputExcel);
    end

    writetable(FrameResolvedTable,outputExcel,'Sheet','FrameResolved');
    writetable(FrameSummaryTable,outputExcel,'Sheet','FrameSummary');
    writetable(ExperimentSettingsTable,outputExcel,'Sheet','ExperimentSettings');
    writematrix(BP_Matrix,outputExcel,'Sheet','BP_Matrix');
    writematrix(MeanIntensity_Matrix,outputExcel,'Sheet','MeanIntensity_Matrix');
    writematrix(STDIntensity_Matrix,outputExcel,'Sheet','STDIntensity_Matrix');
end

if saveMAT

    outputMAT = fullfile(outputDir, ...
        'MultiExperiment_BundlingParameter_AllFrames.mat');

    save(outputMAT, ...
        'BP_Matrix','MeanIntensity_Matrix','STDIntensity_Matrix', ...
        'AverageProfile_ALL','STDProfileAcrossRows_ALL', ...
        'MeanBP','SDBP','SEMBP','NumberValidExperimentsByFrame', ...
        'experimentNumbers','framesToAnalyze','frameInterval_s','time_s', ...
        'sharedROIPosition_px','experimentShift_px', ...
        'ROIPositionUsed_px','ROIBounds_px', ...
        'activationCenters_um','activationRadii_um','pixelSize_um', ...
        'TIFFPaths','NumberAvailableFrames', ...
        'ConditionPassed','ConditionMessage', ...
        'FrameResolvedTable','FrameSummaryTable', ...
        'ExperimentSettingsTable','-v7.3');
end

if saveFigures
    saveVectorAndPNG(figBP,outputDir,'BundlingParameter_AllExperiments');
    saveVectorAndPNG(figHeatmap,outputDir,'BundlingParameter_Heatmap');

    if ~isempty(figDiagnostic)
        saveVectorAndPNG(figDiagnostic,outputDir,'SharedROI_Diagnostic');
    end
end

fprintf('\n============================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('============================================\n');
fprintf('Valid experiments: %d/%d\n',nnz(ConditionPassed),numberExperiments);
fprintf('Output folder:\n%s\n',outputDir);

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function [xStart,yStart,xEnd,yEnd] = rectangleToPixelBounds( ...
        rectanglePosition,numberColumns,numberRows)

    xStart = max(1,round(rectanglePosition(1)));
    yStart = max(1,round(rectanglePosition(2)));

    xEnd = min(numberColumns, ...
        round(rectanglePosition(1)+rectanglePosition(3)));

    yEnd = min(numberRows, ...
        round(rectanglePosition(2)+rectanglePosition(4)));

    if xEnd <= xStart || yEnd <= yStart
        error('The selected rectangle is too small or outside the image.');
    end
end

function displayPurpleImage(axisHandle,imageMatrix,displayPercentiles)

    imagesc(axisHandle,imageMatrix);
    axis(axisHandle,'image');
    set(axisHandle,'YDir','normal');
    colormap(axisHandle,purpleColormap(256));

    finitePixels = imageMatrix(isfinite(imageMatrix));

    if ~isempty(finitePixels)
        limits = prctile(finitePixels,displayPercentiles);
        if limits(2) > limits(1)
            clim(axisHandle,limits);
        end
    end
end

function cmap = purpleColormap(numberColors)

    anchorColors = [ ...
        1.00,1.00,1.00; ...
        0.91,0.84,0.94; ...
        0.66,0.43,0.75; ...
        0.38,0.14,0.49; ...
        0.12,0.03,0.17];

    cmap = interp1( ...
        linspace(0,1,size(anchorColors,1)),anchorColors, ...
        linspace(0,1,numberColors),'pchip');

    cmap = min(max(cmap,0),1);
end

function drawActivationCircles( ...
        axisHandle,centers_px,radii_px,circleColor,lineWidth)

    theta = linspace(0,2*pi,500);

    for circleIndex = 1:numel(radii_px)
        xCircle = centers_px(circleIndex,1) + ...
            radii_px(circleIndex)*cos(theta);
        yCircle = centers_px(circleIndex,2) + ...
            radii_px(circleIndex)*sin(theta);

        plot(axisHandle,xCircle,yCircle,'-', ...
            'Color',circleColor,'LineWidth',lineWidth);
    end
end

function formatAxes(axisHandle,fontSize,blackColor)

    set(axisHandle, ...
        'FontName','Arial','FontSize',fontSize, ...
        'LineWidth',1.5,'TickDir','in', ...
        'TickLength',[0.018,0.018],'Box','on', ...
        'XColor',blackColor,'YColor',blackColor,'Layer','top');

    grid(axisHandle,'off');
end

function saveVectorAndPNG(figureHandle,outputDirectory,baseName)

    if isempty(figureHandle) || ~isgraphics(figureHandle)
        return;
    end

    set(figureHandle,'Renderer','painters');

end
