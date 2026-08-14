%% Paired comparison of stepping fractions in conditions A and B
% Each A-B pair comes from the same independent experiment (n = 3).
clear;
clc;
% close all;

%% Data
% Fraction (stepping emitters/long emitters)
A = [4/22, 5/10, 9/19];
B = [28/106, 21/72, 11/45];
nExperiments = numel(A);

%% Summary statistics
meanA = mean(A);
meanB = mean(B);
sdA = std(A);
sdB = std(B);
pairedDifference = A - B;

%% Paired t-tests
% Two-sided test: A ~= B
[hTwoSided, pTwoSided, ciDifference, statsTwoSided] = ttest(A, B);

% One-sided test: A > B
[hAGreater, pAGreater] = ttest(A, B, 'Tail', 'right');

% One-sided test: B > A
[hBGreater, pBGreater] = ttest(A, B, 'Tail', 'left');

fprintf('Paired measurements (fractions):\n');
fprintf('Experiment %d: A = %.4f, B = %.4f, A-B = %+.4f\n', ...
    [(1:nExperiments); A; B; pairedDifference]);
fprintf('\nA: %.4f +/- %.4f (mean +/- SD)\n', meanA, sdA);
fprintf('B: %.4f +/- %.4f (mean +/- SD)\n', meanB, sdB);
fprintf('Mean paired difference (A-B): %.4f\n', mean(pairedDifference));
fprintf('95%% CI of paired difference: [%.4f, %.4f]\n', ...
    ciDifference(1), ciDifference(2));
fprintf('Two-sided paired t-test: t(%d) = %.4f, p = %.4f, h = %d\n', ...
    statsTwoSided.df, statsTwoSided.tstat, pTwoSided, hTwoSided);
fprintf('One-sided paired t-test (A > B): p = %.4f, h = %d\n', ...
    pAGreater, hAGreater);
fprintf('One-sided paired t-test (B > A): p = %.4f, h = %d\n', ...
    pBGreater, hBGreater);

%% Publication-style paired plot
fig = figure('Color', 'w', 'Position', [100, 100, 520, 520]);
ax = axes(fig);
hold(ax, 'on');

colorA = [0.20, 0.45, 0.75];
colorB = [0.85, 0.33, 0.20];
pairColor = [0.68, 0.68, 0.68];

% Connect measurements from the same experiment.
for experiment = 1:nExperiments
    plot(ax, [1, 2], [A(experiment), B(experiment)], '-', ...
        'Color', pairColor, 'LineWidth', 1.3);
end

% Plot every independent experiment.
scatter(ax, ones(1, nExperiments), A, 70, colorA, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
scatter(ax, 2 * ones(1, nExperiments), B, 70, colorB, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.7);

% Add mean +/- SD. The small horizontal offset keeps these visible beside
% the individual data points.
errorbar(ax, 0.88, meanA, sdA, 'o', 'Color', colorA, ...
    'MarkerFaceColor', 'w', 'MarkerSize', 7, 'LineWidth', 1.8, ...
    'CapSize', 10);
errorbar(ax, 2.12, meanB, sdB, 'o', 'Color', colorB, ...
    'MarkerFaceColor', 'w', 'MarkerSize', 7, 'LineWidth', 1.8, ...
    'CapSize', 10);

% Statistical annotation.
yMax = max([A, B, meanA + sdA, meanB + sdB]);
bracketY = min(0.90, yMax + 0.10);
tickHeight = 0.025;
plot(ax, [1, 1, 2, 2], ...
    [bracketY - tickHeight, bracketY, bracketY, bracketY - tickHeight], ...
    'k-', 'LineWidth', 1.2);
text(ax, 1.5, bracketY + 0.025, sprintf('p = %.3f', pTwoSided), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 11);

ax.XLim = [0.55, 2.45];
ax.YLim = [0, max(1.0, bracketY + 0.12)];
ax.XTick = [1, 2];
ax.XTickLabel = {'No flow', 'With flow'};
ax.YTick = 0:0.2:1;
ax.YTickLabel = compose('%d', 0:20:100);
ax.FontName = 'Arial';
ax.FontSize = 12;
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.Box = 'off';
ylabel(ax, 'stepping / long emitters (%)');

%%