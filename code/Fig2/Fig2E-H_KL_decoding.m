clearvars;
close all;
clc;

%% To test 1: Normalization for the final matrix

%%%
% Normalization Function
function normalized_matrix = normalize_with_threshold(matrix)
    % Fit a normal distribution to the data
    pd = fitdist(matrix(:), 'Normal');
    
    % Extract the mean and standard deviation from the fitted distribution
    mean_val = pd.mu;
    std_dev = pd.sigma;

    % Define threshold values
    upper_threshold = mean_val + 3 * std_dev;
    lower_threshold = mean_val - 3 * std_dev;

    % Replace values larger than upper threshold
    matrix(matrix > upper_threshold) = upper_threshold;

    % Replace values smaller than lower threshold
    matrix(matrix < lower_threshold) = lower_threshold;

    % Normalize the new matrix
    normalized_matrix = (matrix - min(matrix(:))) ./ (max(matrix(:))- min(matrix(:)));
end
%%%

%%%
% Normalization Function
function justnormalize = JustNormalize(matrix)

    % Normalize the new matrix
    justnormalize = (matrix - min(matrix(:))) ./ (max(matrix(:))- min(matrix(:)));
end
%%%

function probMatrix = normalize_matrix(intensityMatrix)
    % Ensure the matrix has non-negative values
    if any(intensityMatrix(:) < 0)
        error('Matrix contains negative values');
    end

    % Calculate the sum of all elements in the matrix
    totalSum = sum(intensityMatrix(:));

    % Normalize the matrix to convert it to a probability distribution
    probMatrix = intensityMatrix / totalSum;
end


function D_KL = calculate_kl_divergence(P, Q)
    % Calculate the KL divergence between two probability matrices P and Q
    % P: First probability matrix
    % Q: Second probability matrix (reference distribution)
    % D_KL: KL divergence matrix

    % Ensure P and Q have the same size
    assert(all(size(P) == size(Q)), 'Matrices P and Q must have the same size');
    
    % Replace zero values in Q to avoid division by zero
    Q(Q == 0) = eps;
    
    % Replace zero values in P to avoid log(0)
    P(P == 0) = eps;
    
    % Calculate the KL divergence
    D_KL = sum(P .* log(P ./ Q), 'all');
end




%% Read setting file

Experiments = 1;
directory = ['../../Experiments/MultiLightSpots_Data&PIV/set' num2str(Experiments) '/'];
%%directory = 'set1/';
files = dir([directory 'PIV_Results/']);
num_files = length(files) - 2;  % Subtract 2 to exclude '.' and '..' directories
start_index = 1;
end_index = num_files;

% find the unit for x and y
txtFiles = dir(fullfile(directory, '*.txt'));
filePath = fullfile(directory, txtFiles(1).name);
fileID = fopen(filePath, 'r');
C = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
lines = C{1,1};
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

% find the light area, light radius
CENTERS  = zeros(12,2);
LRadius = [];
for p = 1:12
    line = lines{end-3-p, 1};
    start_idx = strfind(line, '(');
    end_idx = strfind(line, ')');
    points = [];
    for i = 1:length(start_idx)
        point_str = line(start_idx(i)+1:end_idx(i)-1);
        coords = strsplit(point_str, ',');
        if length(coords) == 2
            point = [str2double(strtrim(coords{1})), str2double(strtrim(coords{2}))];
            points = [points; point];
        end
    end
    center = [0.5*(points(1,1)+points(2,1))*x_unit 0.5*(points(1,2)+points(2,2))*y_unit];
    CENTERS(p,:) = center;
    LRadius(end+1) = (1/(2*sqrt(2)))*sqrt(((points(2,1)-points(1,1))*x_unit)^2+((points(2,2)-points(1,2))*y_unit)^2);
end

X_Centers = CENTERS(:,1);
Y_Centers = CENTERS(:,2);
R_light = mean(LRadius);
R_light  = round(R_light * 1000) / 1000;

%% Initialize the Coord Mesh

tifFiles = dir(fullfile(directory, '*.tif'));
filePath = fullfile(directory, tifFiles(1).name);
I = imread(filePath, 'Index', 1);


unique_x = [1:1:length(I(1,:))]*x_unit;
unique_y = [1:1:length(I(:,1))]*y_unit;
[X, Y] = meshgrid(unique_x, unique_y);

%%% initialize the Intensity matrix for this experiment
Intensity_Array = zeros(length(I(:,1)), length(I(1,:)), end_index, 'uint16');  % Preallocate a 3D array


%% STD Results extracting

for i = start_index:end_index+2
    Intensity_Array(:,:,i) = imread(filePath,'Index',i);
end

Intensity_Std = std(double(Intensity_Array), 0, 3);


UnINTENSITY_Std = Intensity_Std;
Intensity_Std = normalize_with_threshold(Intensity_Std);


%% Experiment Result extracting

Row_coords = [53.9 73 92.45 111.9 133 152.25 169.875 187.5];
Colum_coords = [65 87 108 128];
RUNNING = 0;
Intensity_Std_Average = [];

Inward_Index = [5 6 7 8 13 17 19 23 24 29 30 31 32];

for m = 1:length(Row_coords)
    for n = 1:length(Colum_coords)
    RUNNING = RUNNING+1;
    %if ismember(RUNNING,Inward_Index) == 1
    Center_Point = [Colum_coords(n), Row_coords(m)];
    center_x = Center_Point(1);
    center_y = Center_Point(2);
        
    inner_radius = R_light*0;
    outer_radius = R_light*1;

    distance_from_center = sqrt((X - center_x).^2 + (Y - center_y).^2);
    annular_mask = (distance_from_center >= inner_radius) & (distance_from_center <= outer_radius);
    
    % Extract the values within the annulus using the mask
    
    values_in_annulus = Intensity_Std(annular_mask);

    Intensity_Std_Average(end+1) = mean(values_in_annulus);
    %else
    %Intensity_Std_Average(end+1) = 0;
    %end
    end
end

INTENSITY_STD = Intensity_Std_Average;
INTENSITY_STD = reshape(INTENSITY_STD,4,8);
INTENSITY_STD = INTENSITY_STD';

INTENSITY_VARIANCE = INTENSITY_STD;
INTENSITY_VARIANCE = JustNormalize(INTENSITY_VARIANCE);
INTENSITY_VARIANCE = normalize_matrix(INTENSITY_VARIANCE);


 
%% All Possibilities

Inward_Index = [5 6 7 8 13 17 19 23 24 29 30 31 32];

n = length(Inward_Index); % How many positions
Number_LightSpots = 1:1:n; % How many light spots are activated

ConfigurationAndSSIM = {};

filename = 'StdConfigurations.txt';
% Open the text file in append mode
fileID = fopen(filename, 'a');


for k = Number_LightSpots % How many light points are activated
index_combinations = nchoosek(1:n, k);
    for o = 1:length(index_combinations(:,1))
        Possible_STD = zeros(1,32);
 
        disp([num2str(k) 'Light Spots activated___' num2str(o) 'th Configuration'])
  
        for w = 1:k
            Possible_STD(Inward_Index(index_combinations(o,w))) = 1;
        end
        Possible_STD = reshape(Possible_STD,4,8);
        Possible_STD = Possible_STD';
        ConfigurationAndSSIM{end+1,1} = Inward_Index(index_combinations(o,:));
        Possible_STD = normalize_matrix(Possible_STD);
        
        ConfigurationAndSSIM{end,2} = calculate_kl_divergence(Possible_STD,INTENSITY_VARIANCE);

        % Prepare the data to write
        dataToWrite = {ConfigurationAndSSIM{end,1},ConfigurationAndSSIM{end,2}};
        % Format the data and write it to the file
        fprintf(fileID, '%s\t%f\n', dataToWrite{:});
    end
end
fclose(fileID);

%% Answer

answer = [5 6 7 8 13 17 19 23 24 29 30 32];
ANS_ConfigurationAndSSIM = {};

Possible_STD = zeros(1,32);


for i = answer
    Possible_STD(i) = 1;
end

Possible_STD_ANS = Possible_STD;
Possible_STD_ANS = reshape(Possible_STD_ANS,4,8);
Possible_STD_ANS = Possible_STD_ANS';
Possible_STD_ANS = normalize_matrix(Possible_STD_ANS);

ANS_ConfigurationAndSSIM{end+1,1} = answer;
ANS_ConfigurationAndSSIM{end,2} = calculate_kl_divergence(INTENSITY_VARIANCE,Possible_STD_ANS);


%% Pick the best configuration
% Define the experimental data
% Each row corresponds to an experiment, and each column corresponds to a parameter
experiments = cell2mat(ConfigurationAndSSIM(:, 2));

% Number of experiments
num_experiments = size(experiments, 1);

best_experiment_index = find(experiments==min(experiments));

% Get the best experiment
best_SSIM = experiments(best_experiment_index, :)
best_configuration = cell2mat(ConfigurationAndSSIM(best_experiment_index,1))

%% Plot PIV

left_x = 65 - 3 * R_light;
right_x = 128 + 3 * R_light;
upper_y = 187.5 + 3 * R_light;
lower_y = 53.9 - 3 * R_light;

[~, x_start_idx] = min(abs(unique_x - left_x));
[~, y_start_idx] = min(abs(unique_y - lower_y));
[~, x_end_idx] = min(abs(unique_x - right_x));
[~, y_end_idx] = min(abs(unique_y - upper_y));


Intensity_Std = Intensity_Std(y_start_idx:y_end_idx,x_start_idx:x_end_idx);
UnINTENSITY_Std = UnINTENSITY_Std(y_start_idx:y_end_idx,x_start_idx:x_end_idx);


unique_x = unique_x(x_start_idx:x_end_idx);
unique_y = unique_y(y_start_idx:y_end_idx);

X = X(y_start_idx:y_end_idx,x_start_idx:x_end_idx);
Y = Y(y_start_idx:y_end_idx,x_start_idx:x_end_idx);

%%
Row_coords = [53.9 73 92.45 111.9 133 152.25 169.875 187.5];
Colum_coords = [65 87 108 128];
Center_Point = {};
for i = 1:length(Row_coords)
    for j = 1:length(Colum_coords)
        Center_Point{end+1} = [Colum_coords(j),max(unique_y)+min(unique_y)-Row_coords(i)];
    end
end


%%
FZ = 25;
figure(1)
h = pcolor(X, max(unique_y)+min(unique_y)-Y, Intensity_Std);
hold on
set(h, 'EdgeColor', 'none'); 
shading interp
cb = colorbar;
cb.Location = 'westoutside';
cb.Label.FontSize = FZ;
pos = get(cb, 'Position');
pos(1) = pos(1)+0.132; % Adjust the value to move the colorbar closer to the plot
set(cb, 'Position', pos);
cb.Label.String = 'Normalized STD';
axis equal
xlim([min(unique_x) max(unique_x)])
ylim([min(unique_y) max(unique_y)])
%xlabel('X (µm)','FontSize',FZ)
%ylabel('Y (µm)','FontSize',FZ)
%title('Normalized STD')
axis off
clim([0 1])

RUNNING = 0;
for i = Row_coords
    for j = Colum_coords
        RUNNING = RUNNING+1;
    center = Center_Point{RUNNING};
    theta = linspace(0, 2*pi, 100);
    radius = R_light;
    X_Circle = center(1) + radius * cos(theta);
    Y_Circle = center(2) + radius * sin(theta);
    plot(X_Circle,Y_Circle,'r-',LineWidth=2.5)
    end
end

hold off
set(gca,'fontsize',FZ,'linewidth',1,'ticklength',[0.02 0.02])
%exportgraphics(gcf,'STDField.eps','contenttype','vector')
%saveas(gcf,"STDField.png")
%%
Row_coords = [53.9 73 92.45 111.9 133 152.25 169.875 187.5];
Colum_coords = [65 87 108 128];
Center_Point = {};
for i = 1:length(Row_coords)
    for j = 1:length(Colum_coords)
        Center_Point{end+1} = [Colum_coords(j),max(unique_y)+min(unique_y)-Row_coords(i)];
    end
end

FZ = 15;
figure(4)
h = pcolor(X, max(unique_y)+min(unique_y)-Y, UnINTENSITY_Std);
hold on
set(h, 'EdgeColor', 'none'); 
shading interp
cb = colorbar;
cb.Label.String = 'Intensity STD';
axis equal
xlim([min(unique_x) max(unique_x)])
ylim([min(unique_y) max(unique_y)])
xlabel('X (µm)','FontSize',FZ)
ylabel('Y (µm)','FontSize',FZ)
title('Raw STD')

RUNNING = 0;
for i = Row_coords
    for j = Colum_coords
        RUNNING = RUNNING+1;
    center = Center_Point{RUNNING};
    theta = linspace(0, 2*pi, 100);
    radius = R_light;
    X_Circle = center(1) + radius * cos(theta);
    Y_Circle = center(2) + radius * sin(theta);
    plot(X_Circle,Y_Circle,'r-',LineWidth=1.2)
    end
end

hold off
set(gca,'fontsize',FZ,'linewidth',1,'ticklength',[0.02 0.02])
%%
% Plot Answer
figure(2)
h = pcolor(X, max(unique_y)+min(unique_y)-Y, Intensity_Std);
hold on
set(h, 'EdgeColor', 'none'); 
shading interp
cb = colorbar;
cb.Label.String = 'Normalized Intensity STD';
axis equal
xlim([min(unique_x) max(unique_x)])
ylim([min(unique_y) max(unique_y)])
xlabel('X (µm)','FontSize',FZ)
ylabel('Y (µm)','FontSize',FZ)
title('Answer',FontSize=20)

Row_coords = [53.9 73 92.45 111.9 133 152.25 169.875 187.5];
Colum_coords = [65 87 108 128];
for i = answer
    center = Center_Point{i};
    theta = linspace(0, 2*pi, 100);
    radius = R_light;
    X_Circle = center(1) + radius * cos(theta);
    Y_Circle = center(2) + radius * sin(theta);
    plot(X_Circle,Y_Circle,'r-',LineWidth=1.2)
end
hold off
set(gca,'fontsize',FZ,'linewidth',1,'ticklength',[0.02 0.02])
%%

% Plot the best case
figure(3)
h = pcolor(X, max(unique_y)+min(unique_y)-Y, Intensity_Std);
hold on
set(h, 'EdgeColor', 'none'); 
shading interp
cb = colorbar;
cb.Label.String = 'Normalized Intensity STD';
axis equal

xlim([min(unique_x) max(unique_x)])
ylim([min(unique_y) max(unique_y)])
xlabel('X (µm)','FontSize',FZ)
ylabel('Y (µm)', 'FontSize',FZ)

title('Normalized STD')


configuration = best_configuration;
for i = configuration
    center = Center_Point{i};
    theta = linspace(0, 2*pi, 100);
    radius = R_light;
    X_Circle = center(1) + radius * cos(theta);
    Y_Circle = center(2) + radius * sin(theta);
    plot(X_Circle,Y_Circle,'r-',LineWidth=2)
end
hold off
set(gca,'fontsize',FZ,'linewidth',1,'ticklength',[0.02 0.02])





%%
figure(5)
hold on

axis equal
xlim([min(unique_x) max(unique_x)])
ylim([min(unique_y) max(unique_y)])
%xlabel('X (µm)','FontSize',FZ)
%ylabel('Y (µm)','FontSize',FZ)
title('best')
axis off

Row_coords = [53.9 73 92.45 111.9 133 152.25 169.875 187.5];
Colum_coords = [65 87 108 128];
for i = 1:32
    center = Center_Point{i};
    theta = linspace(0, 2*pi, 100);
    radius = R_light;
    X_Circle = center(1) + radius * cos(theta);
    Y_Circle = center(2) + radius * sin(theta);
    if ismember(i,best_configuration) == 1
        fill(X_Circle,Y_Circle,'r-',LineWidth=1.2)
    else
        fill(X_Circle,Y_Circle,[0.7 0.7 0.7],LineWidth=1.2)
    end
end

box on
hold off
set(gca,'fontsize',FZ,'linewidth',1,'ticklength',[0.02 0.02])



%%
figure(6)
hold on

axis equal
xlim([min(unique_x) max(unique_x)])
ylim([min(unique_y) max(unique_y)])
%xlabel('X (µm)','FontSize',FZ)
%ylabel('Y (µm)','FontSize',FZ)
%title('Another random configuration')
axis off

Row_coords = [53.9 73 92.45 111.9 133 152.25 169.875 187.5];
Colum_coords = [65 87 108 128];
for i = 1:32
    center = Center_Point{i};
    theta = linspace(0, 2*pi, 100);
    radius = R_light;
    X_Circle = center(1) + radius * cos(theta);
    Y_Circle = center(2) + radius * sin(theta);
    if ismember(i,ConfigurationAndSSIM{8102,1}) == 1
        fill(X_Circle,Y_Circle,'y-',LineWidth=1.2)
    else
        fill(X_Circle,Y_Circle,[0.7 0.7 0.7],LineWidth=1.2)
    end
end

box off
hold off
set(gca,'fontsize',FZ,'linewidth',1,'ticklength',[0.02 0.02])

exportgraphics(gcf, 'competing.pdf', 'ContentType', 'vector');

%%

%% Find the smallest and second-smallest KL configurations

% Column 2 contains the KL-divergence values
KL_values = cell2mat(ConfigurationAndSSIM(:,2));

% Sort KL values in ascending order
[sorted_KL, sorted_idx] = sort(KL_values, 'ascend');

% Smallest KL result
smallest_id = sorted_idx(1);
smallest_KL = sorted_KL(1);
smallest_configuration = ConfigurationAndSSIM{smallest_id,1};

% Second-smallest KL result
second_smallest_id = sorted_idx(2);
second_smallest_KL = sorted_KL(2);
second_smallest_configuration = ConfigurationAndSSIM{second_smallest_id,1};

% Display results
fprintf('\n===== Best two configurations =====\n');

fprintf('Smallest KL:\n');
fprintf('  ID = %d\n', smallest_id);
fprintf('  KL = %.10g\n', smallest_KL);
fprintf('  Configuration = ');
fprintf('%d ', smallest_configuration);
fprintf('\n');

fprintf('Second-smallest KL:\n');
fprintf('  ID = %d\n', second_smallest_id);
fprintf('  KL = %.10g\n', second_smallest_KL);
fprintf('  Configuration = ');
fprintf('%d ', second_smallest_configuration);
fprintf('\n');