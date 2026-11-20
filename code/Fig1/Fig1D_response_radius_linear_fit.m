clearvars;
close all;
clc;

filename = 'SortedData.xlsx'; 
data = readmatrix(filename);

% Extract columns
x = data(:, 2);  % Assuming the first column is x
y = data(:, 3);  % Assuming the second column is y
weights = data(:, 4);  % Assuming the third column is weights
%weights = ones(28,1);

% Perform linear regression with weights
mdl = fitlm(x, y, 'Weights', weights);

% Display summary of the linear model
disp(mdl);

% Extract fitting model
FittingModel = mdl.Coefficients;

% filename = ['FittingResults.xlsx'];
% writetable(FittingModel, filename);
% 
% LRNo = [154.5, 126, 117.5, 122]*0.1;
% RRNo = [299, 318, 370, 331]*0.1;

%% Plot

rgb_Purple = [135, 85, 155]/255;
rgb_Yellow = [227,207,169]/255;

FZ = 25;
LZ = 2;
% Plot scatter plot
figure(1);
scatter(x, y,'MarkerEdgeColor',rgb_Purple,'SizeData',150,LineWidth=3, Marker='diamond')
hold on;
% Plot regression line
x_range = 0:0.1:40;
y_pred = predict(mdl, x_range');
plot(x_range, y_pred, 'k--', 'LineWidth', 5);
%scatter(LRNo, RRNo, 'ro','LineWidth',3, 'SizeData', 80)
hold off;
xlim([0 40])
ylim([15 100])
box on
axis square
xlabel('Laser Radius (µm)','FontSize',FZ)
ylabel('Response Radius by PIV (µm)','FontSize',FZ)
%leg1=legend('Raw data', 'Fitting Model');
%set(leg1,'box','off','location','northwest','fontsize',15,'Numcolumn',1)
%title('Relationship between Light And Activation')
set(gca,'fontsize',FZ,'linewidth',1,'ticklength',[0.02 0.02])
saveas(1,'FittingSingleLightSpot.png')
exportgraphics(gcf,'e2.eps','contenttype','vector')