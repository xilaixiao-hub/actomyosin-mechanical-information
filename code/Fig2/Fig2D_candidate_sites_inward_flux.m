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
    upper_threshold = mean_val + 100 * std_dev;
    lower_threshold = mean_val - 100 * std_dev;

    % Replace values larger than upper threshold
    matrix(matrix > upper_threshold) = upper_threshold;

    % Replace values smaller than lower threshold
    matrix(matrix < lower_threshold) = lower_threshold;

    % Normalize the new matrix
    normalized_matrix = (matrix - min(matrix(:))) ./ (max(matrix(:))- min(matrix(:)));
end
%%%

EXPERIMENTS = 1:1;

for EXPERIMENT = EXPERIMENTS

    
%% Read setting file

Experiments = 1;
directory = ['../../Experiments/MultiLightSpots_Data&PIV/set' num2str(Experiments) '/'];
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

%% Initialize the Coord Mesh

filename = fullfile([directory 'PIV_Results'], 'PIV_t1-t3.txt');
data = load(filename);
x = x_unit*data(:, 1);

y = y_unit*data(:, 2);

unique_x = unique(x);
unique_y = unique(y);
[X, Y] = meshgrid(unique_x, unique_y);

DISPLACEMENT_X = zeros(length(unique_y), length(unique_x));
DISPLACEMENT_Y = zeros(length(unique_y), length(unique_x));

%% PIV Results extracting
net_flow = zeros(32,1);

for i = start_index:end_index
    filename = fullfile([directory 'PIV_Results'], ['PIV_t' num2str(i) '-t' num2str(i+2) '.txt']);
    x = x_unit*data(:, 1);
    y = y_unit*data(:, 2);

    if exist(filename, 'file') == 2
    data = load(filename);
    u = x_unit*data(:, 13)/2;      % (After interpolated, ux0)
    v = y_unit*data(:, 14)/2;      % (After interpolated, uy0)
    else
    u = zeros(length(unique_x)*length(unique_y),1);
    v = zeros(length(unique_x)*length(unique_y),1);
    end

    % Creating grids for X, Y for Magnitude Plot
    U = griddata(x, y, u, X, Y);
    V = griddata(x, y, v, X, Y);

%% To test 2: Median Replacement for each PIV result
    %%%     
    % Median Replacement for U
    matrix = U;
    % Calculate the median of the entire matrix
    median_value = median(matrix(:));
    
    % Find indices of extreme values (outliers)
    outlier_indices = abs(matrix - median_value) > 2 * std(matrix(:)); % Adjust the threshold multiplier as needed
    
    % Interpolate extreme values using median of surrounding values
    for y = 1:numel(matrix)
        if outlier_indices(y)
            [row, col] = ind2sub(size(matrix), y);
            
            % Find surrounding values
            surrounding_values = matrix(max(row-1,1):min(row+1,end), max(col-1,1):min(col+1,end));
            
            % Calculate median of surrounding values
            interpolated_value = median(surrounding_values(:));
            
            % Replace extreme value with interpolated value
            matrix(row, col) = interpolated_value;
        end
    end
    U = matrix;
    %%%

    %%%
    % Median Replacement for V
    matrix = V;
    % Calculate the median of the entire matrix
    median_value = median(matrix(:));
    
    % Find indices of extreme values (outliers)
    outlier_indices = abs(matrix - median_value) > 2 * std(matrix(:)); % Adjust the threshold multiplier as needed
    
    % Interpolate extreme values using median of surrounding values
    for y = 1:numel(matrix)
        if outlier_indices(y)
            [row, col] = ind2sub(size(matrix), y);
            
            % Find surrounding values
            surrounding_values = matrix(max(row-1,1):min(row+1,end), max(col-1,1):min(col+1,end));
            
            % Calculate median of surrounding values
            interpolated_value = median(surrounding_values(:));
            
            % Replace extreme value with interpolated value
            matrix(row, col) = interpolated_value;
        end
    end
    V = matrix;
    %%%


Row_coords = [53.9 73 92.45 111.9 133 152.25 169.875 187.5];
Colum_coords = [65 87 108 128];
RUNNING = 0;
InwardFlow = {};


for m = 1:length(Row_coords)
    for n = 1:length(Colum_coords)
    RUNNING = RUNNING+1
    Center_Point = [Colum_coords(n), Row_coords(m)];
    center_x = Center_Point(1);
    center_y = Center_Point(2);        
    
    
        % Define the circle parameters
    x_c = center_x; % x-coordinate of the center
    y_c = center_y; % y-coordinate of the center
    r = R_light;   % radius of the circle
    
    % Discretize the circle
    num_points = 100; % number of points on the circle
    theta = linspace(0, 2*pi, num_points);
    x = x_c + r * cos(theta);
    y = y_c + r * sin(theta);
    
    % Load or define the velocity field
    % Assuming velocity field is given as u(x, y) and v(x, y)
    % Example velocity field: u = -y, v = x (a rotational field)
    u = U;
    v = V;

    % Calculate the net flow
    for i = 1:num_points-1
    % Midpoint of the segment
    x_mid = (x(i) + x(i+1)) / 2;
    y_mid = (y(i) + y(i+1)) / 2;
    
    % Normal vector at the midpoint
    normal_vector = [x_mid - x_c, y_mid - y_c];
    normal_vector = normal_vector / norm(normal_vector);
    
    % Velocity at the midpoint
    u_mid = interp2(X, Y, U, x_mid, y_mid, 'linear');
    v_mid = interp2(X, Y, V, x_mid, y_mid, 'linear');

    % Normal component of the velocity
    v_normal = u_mid * normal_vector(1) + v_mid * normal_vector(2);
    
    % Length of the segment
    segment_length = sqrt((x(i+1) - x(i))^2 + (y(i+1) - y(i))^2);
    
    % Contribution to the net flow
    net_flow(RUNNING) = net_flow(RUNNING) + v_normal * segment_length;
    end

    % All_Radius = 0.1:0.1:R_light;
    % theta = linspace(0, 2*pi, 100);
    % InwardFlow{end+1,1} = [];
    % 
    % for radius = All_Radius
    % % Generate points on the circumference of the circle
    % X_Circle = center_x + radius * cos(theta);
    % Y_Circle = center_y + radius * sin(theta);
    % 
    % % Assuming you have your data in the form of arrays x_data, y_data, and displacement_data
    % % Use griddata to interpolate displacement values at the points on the circle
    % Displacement_X_circle = griddata(unique_x, unique_y, DISPLACEMENT_X, X_Circle, Y_Circle);
    % Displacement_Y_circle = griddata(unique_x, unique_y, DISPLACEMENT_Y, X_Circle, Y_Circle);
    % % Average the interpolated displacement values
    % Displacement_Radial_circle = Displacement_X_circle.*cos(theta)+Displacement_Y_circle.*sin(theta);
    % InwardFlow{end,1}(end+1) = -mean(Displacement_Radial_circle);
    % end
    % 
    % figure;
    % plot(All_Radius,InwardFlow{end,1})
    % title([num2str(RUNNING)])
    end
end


    DISPLACEMENT_X = DISPLACEMENT_X + U;
    DISPLACEMENT_Y = DISPLACEMENT_Y + V;
end
end
%%
DISPLACEMENT = sqrt(DISPLACEMENT_X.^2+DISPLACEMENT_Y.^2);
DISPLACEMENT = normalize_with_threshold(DISPLACEMENT);

left_x = 65 - 3 * R_light;
right_x = 128 + 3 * R_light;
upper_y = 187.5 + 3 * R_light;
lower_y = 53.9 - 3 * R_light;

[~, x_start_idx] = min(abs(unique_x - left_x));
[~, y_start_idx] = min(abs(unique_y - lower_y));
[~, x_end_idx] = min(abs(unique_x - right_x));
[~, y_end_idx] = min(abs(unique_y - upper_y));

% PIV results: DISPLACEMENT field after normalization


DISPLACEMENT = DISPLACEMENT(y_start_idx:y_end_idx,x_start_idx:x_end_idx);
DISPLACEMENT_X = DISPLACEMENT_X(y_start_idx:y_end_idx,x_start_idx:x_end_idx);
DISPLACEMENT_Y = DISPLACEMENT_Y(y_start_idx:y_end_idx,x_start_idx:x_end_idx);

unique_x = unique_x(x_start_idx:x_end_idx);
unique_y = unique_y(y_start_idx:y_end_idx);

X = X(y_start_idx:y_end_idx,x_start_idx:x_end_idx);
Y = Y(y_start_idx:y_end_idx,x_start_idx:x_end_idx);
answer = [1 2 3 4 9 13 15 19 20 25 26 28]+4;

INDEX = find(net_flow < 0);
INDEX = INDEX';


%%
FZ = 25;
LZ = 1.2;
rgb_Green = [92, 173, 85]/255;           % Typical 0–255 RGB (a reddish-orange)
rgb_Blue = [75, 167, 227]/255;
rgb_Red = [150, 26, 73]/255;
rgb_Yellow = [227,207,169]/255;
figure(1)

colorscheme=othercolor('PuOr11',256);
colorscheme = colorscheme(129:256,:);
colormap((colorscheme));

h = pcolor(X, max(unique_y)+min(unique_y)-Y, DISPLACEMENT);
hold on
set(h, 'EdgeColor', 'none');

cb = colorbar;
cb.Location = 'westoutside';
cb.Label.FontSize = FZ;
pos = get(cb, 'Position');
pos(1) = pos(1)+0.132; % Adjust the value to move the colorbar closer to the plot
set(cb, 'Position', pos);
cb.Label.String = 'Normalized Displacement';

Row_coords = [53.9 73 92.45 111.9 133 152.25 169.875 187.5];
Colum_coords = [65 87 108 128];
Center_Point = {};
for i = 1:length(Row_coords)
    for j = 1:length(Colum_coords)
        Center_Point{end+1} = [Colum_coords(j),max(unique_y)+min(unique_y)-Row_coords(i)];
    end
end

for i = INDEX
    center = Center_Point{i};
    theta = linspace(0, 2*pi, 100);
    radius = R_light*0.8;

    X_Circle = center(1) + radius * cos(theta);
    Y_Circle = center(2) + radius * sin(theta);

    plot(X_Circle,Y_Circle,'k-',LineWidth=3)
end

for i = answer
    center = Center_Point{i};
    theta = linspace(0, 2*pi, 100);
    radius = R_light*1.2;

    X_Circle = center(1) + radius * cos(theta);
    Y_Circle = center(2) + radius * sin(theta);

    plot(X_Circle,Y_Circle,'-','Color',rgb_Red,LineWidth=3)
end

spacing_factor = 1.2;
spaced_x = unique_x(1):spacing_factor*(unique_x(end)-unique_x(1))/(length(unique_x)-1):unique_x(end);
spaced_y = unique_y(1):spacing_factor*(unique_y(end)-unique_y(1))/(length(unique_y)-1):unique_y(end);
[Spaced_X, Spaced_Y] = meshgrid(spaced_x, spaced_y);
Spaced_U = griddata(unique_x, unique_y, DISPLACEMENT_X, Spaced_X, Spaced_Y);
Spaced_V = griddata(unique_x, unique_y, DISPLACEMENT_Y, Spaced_X, Spaced_Y);
% Create the quiver plot with downsampled and averaged data
quiver(Spaced_X, Spaced_Y, Spaced_U, Spaced_V, 'Color', 'black', 'LineWidth', 1,'ShowArrowHead', 'on', 'AutoScaleFactor', 2, 'AutoScale', 'on');

hold off
%xlabel('X (µm)','FontSize',FZ)
%ylabel('Y (µm)','FontSize',FZ)
axis off
axis equal
xlim([min(unique_x) max(unique_x)]);
ylim([min(unique_y) max(unique_y)]);
%title('Regions with inward flow')
set(gca,'fontsize',FZ,'linewidth',1,'ticklength',[0.02 0.02])
saveas(gcf, "InwardFlow.png")