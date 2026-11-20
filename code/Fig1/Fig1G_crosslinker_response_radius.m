%% Two-sample T-test and Boxplot Visualization
clc; clear; close all;

% Define font size
FZ = 25; 

% Example data
% From0rTo05r = [
% 0.77181172
% 0.803962658
% 0.768515055
% 0.753180232
% 0.703098045
% 0.579150085
% 0.827837576
% 0.760782929
% ];
% From0rTo05r = [
% 0.77181172
% 0.803962658
% 0.805775069
% 0.703098045
% 0.768515055
% 0.753180232
% 0.756529511
% 0.778221434
% 0.634795321
% 0.688491987
% 0.579150085
% 0.713032364
% 0.827837576
% 0.760782929
% ]
% From05rTo10r = [
% 0.772228422
% 0.601936867
% 0.704472919
% 0.709802566
% 0.659686984
% 0.789858123
% 0.543368854
% 0.698316389
% 0.578141729
% ];

% Data
From05rTo10r = [0.519244709
0.593165955
0.544284957
0.464126532
0.461088442
0.541749073
0.525432261
0.644381925
0.500884561
0.359604205
0.398974696
0.394973132
0.503777852
0.660278656
0.742846328
0.584637089
0.557435711
0.45763181
0.555262164
0.491286373];
From05rTo10r = From05rTo10r';

From0rTo05r = [0.84312122
0.520407063
0.550488848
0.442284276
0.226227813
1.086810175
0.692041816
0.41610271
1.353250121
0.388399818
1.255152384
0.350694612
0.522197513
0.742333287
0.952708867
0.843911411
1.094647844
0.648060984
0.708818713
0.532159135
0.919498587
0.498187183];
From0rTo05r = From0rTo05r';


% From0rTo05r = [34,29,32,31,32,36,31,36,43,30,28];
% %From0rTo05r = From0rTo05r';
% 
% From05rTo10r = [20,35,51,42,73,42,46,61,61];
% %From05rTo10r = From05rTo10r';

% Combine all data into one vector
allData = [From0rTo05r, From05rTo10r];

% Create group labels. Using strings is good practice when specifying positions.
group = [repmat({'0'}, length(From0rTo05r), 1); ...
         repmat({'1/50'}, length(From05rTo10r), 1)];

% Create a figure for the plot
figure;
hold on;

% --- Define the positions for the plots ---
boxPositions = [0.7, 1.4];

% Create a boxplot for the data using the specified positions
boxplot(allData, group,'Positions', boxPositions, 'Colors', 'k', 'Widths', 0.5, 'Whisker', 1.5, 'Symbol', '');

% Set transparency level
alpha_value = 0.5;
rgb_Purple = [135, 85, 155]/255;

% --- Plot swarm charts using the SAME positions as the boxplot ---
% This ensures the swarm points are centered on the boxes.
swarmchart(boxPositions(1) * ones(size(From0rTo05r)), From0rTo05r, 60, rgb_Purple, 'filled', 'MarkerFaceAlpha', alpha_value);
swarmchart(boxPositions(2) * ones(size(From05rTo10r)), From05rTo10r, 60, 'k', 'filled', 'MarkerFaceAlpha', alpha_value);

% Labels and formatting
% Set the x-axis ticks to the box positions and apply labels
set(gca, 'XTick', boxPositions, 'XTickLabel', {'0', '1/50'}, 'FontSize', 25);
ylabel('R', 'FontSize', 25);
ylabel('$R_{rp}\ (\mu m)$','Interpreter','latex')
ylim([0 1.5]);

% Add grid and box for better visualization
box on;
axis square;       % makes width = height in physical units
pbaspect([0.4 1 1]); % same thing explicitly

hold off;


% Perform t-test
[h, p] = ttest2(From0rTo05r, From05rTo10r);

% Title with p-value
%title(['Two-Sample t-test: p = ', num2str(p, '%.3f')], 'FontSize', 20);

% Export the figure
exportgraphics(gcf, 'CV.eps', 'ContentType', 'vector');
