%% compare_two_samples_box_cloud.m
% Box + jittered "cloud" (no broken y-axis), with stats + vector export.
% Matches the style we used earlier: black vs purple, transparent dots,
% thick box lines, square-ish aspect, export as vector EPS.

clear; clc; close all;

%% ================== 1) Input data ==================
% --- Pick ONE pair (comment out the other pair) ---

% Pair A
data1 = [1.272998165, 1.285628441, 1.254042543, 1.221386122, 1.238230326, ...
         1.265482711, 1.287614901, 1.295667606, 1.063680862, 1.139551824, 1.16440598];
data2 = [1.685073104, 1.639798527, 1.594017025, 1.474259065, 1.410558665, ...
         2.700623444, 1.395828472, 2.774614177, 1.388816007, 1.5252, 1.4352];
%%
data1 = [1.25821623
1.269725021
1.097672504
1.139640781
1.348093238
1.273806736
1.246297684
1.27361233
1.315990512
1.221807389
1.309598978
1.221804926
1.25991771
1.115873417
1.064844109
1.436694096
1.401681375
1.528223975];

data2 = [2.285903603
1.783187957
1.862308056
1.784813113
1.802463801
1.647915532
1.502730849
2.603713586
2.040836792
1.638450911
1.609077567
1.548170878
1.562801353
1.30790362
1.60737576
2.798386036
3.923477794
3.215720112];
%%

% % Pair B (uncomment to use)
% data1 = [0.84312122 0.520407063 0.550488848 0.442284276 0.226244329 1.086810175 ...
%          0.692041816 0.41610271 1.353250121 0.388399818 1.255152384 0.350694612 ...
%          0.522197513 0.742333287 0.952708867 0.843911411 1.094647844 0.648060984 ...
%          0.708818713 0.532159135 0.919498587 0.498187183];
% data2 = [0.519244709 0.593165955 0.544284957 0.464126532 0.461088442 0.541749073 ...
%          0.525432261 0.644381925 0.500884561 0.359604205 0.398974696 0.394973132 ...
%          0.503777852 0.660278656 0.742846328 0.584637089 0.557435711];

groupNames = {'0','1/50'};

%% ================== 2) Stats ==================
% Welch t-test (safer if variances differ)
[~, p_value] = ttest2(data1, data2, 'Vartype','unequal');

m1 = mean(data1); m2 = mean(data2);
s1 = std(data1,0); s2 = std(data2,0);
d_prime = abs(m2 - m1) / sqrt(0.5*(s1^2 + s2^2));

fprintf('n1=%d  n2=%d\n', numel(data1), numel(data2));
fprintf('Mean1=%.5f  Mean2=%.5f\n', m1, m2);
fprintf('Std1 =%.5f  Std2 =%.5f\n', s1, s2);
fprintf('Welch ttest2 p=%.4g\n', p_value);
fprintf('d''=%.3f\n', d_prime);

%% ================== 3) Plot knobs (style) ==================
alphaPts    = 0.55;   % dot transparency
dotSize     = 90;     % dot size
jitterWidth = 0.18;   % horizontal spread

FZ = 26;              % axis font size
LW = 2.3;             % boxplot line width

rgb_Purple = [135, 85, 155]/255;

% y-limits control
USE_CUSTOM_YLIM = false;
customYlim = [0 1.5];

% auto y-pad as fraction of range
ylimPadFrac = 0.08;

% output name
outBase = 'box_cloud_no_ybreak';

%% ================== 4) Build figure ==================
fig = figure('Color','w'); hold on;

allData = [data1(:); data2(:)];
group   = [ones(numel(data1),1); 2*ones(numel(data2),1)];

% Boxplot (black outlines, no outlier symbols)
boxplot(allData, group, ...
    'Colors','k', ...
    'Widths',0.50, ...
    'Whisker',1.5, ...
    'Symbol','');

% Thicken boxplot lines (some versions also catch median/whiskers)
set(findobj(gca,'Type','Line'), 'LineWidth', LW);

% Jittered scatter (cloud)
rng(0); % reproducible jitter
x1 = 1 + (rand(numel(data1),1)-0.5)*2*jitterWidth;
x2 = 2 + (rand(numel(data2),1)-0.5)*2*jitterWidth;

scatter(x1, data1, dotSize, 'filled', ...
    'MarkerFaceColor', 'k', ...
    'MarkerFaceAlpha', alphaPts, ...
    'MarkerEdgeAlpha', 0);

scatter(x2, data2, dotSize, 'filled', ...
    'MarkerFaceColor', rgb_Purple, ...
    'MarkerFaceAlpha', alphaPts, ...
    'MarkerEdgeAlpha', 0);

% Axes / labels
set(gca, 'XTick', [1 2], 'XTickLabel', groupNames, ...
    'FontSize', FZ, 'LineWidth', 1.2, 'TickLength', [0.02 0.02]);

ylabel('STD in laser region', 'FontSize', FZ);

% Square-ish look (adjustable)
axis square;
pbaspect([0.4 1 1]);  % change first number to make it wider (e.g., [1.3 1 1])

box on;

% y-limits
if USE_CUSTOM_YLIM
    ylim(customYlim);
else
    dmin = min(allData(isfinite(allData)));
    dmax = max(allData(isfinite(allData)));
    dr = max(dmax - dmin, eps);
    ylim([dmin - ylimPadFrac*dr, dmax + ylimPadFrac*dr]);
end

% Optional annotation (uncomment if you want it on-plot)
% text(0.05, 0.92, sprintf("d' = %.2f", d_prime), 'Units','normalized', ...
%     'FontSize', 20, 'Color','k');
% text(0.05, 0.84, sprintf("p = %.3g", p_value), 'Units','normalized', ...
%     'FontSize', 18, 'Color','k');

%% ================== 5) Export ==================
set(gcf, 'Renderer', 'painters');  % best for EPS vector
drawnow;

exportgraphics(gcf, [outBase '.eps'], 'ContentType','vector');
exportgraphics(gcf, [outBase '.png'], 'Resolution', 300);

fprintf('Saved: %s.eps (vector) and %s.png\n', outBase, outBase);