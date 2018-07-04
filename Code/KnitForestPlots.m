addpath('altmany-export_fig-04ca93c')

close all

fontSize = 8.5;

aa = imread('../Figures/AA_p_Forest.tiff');
bb = imread('../Figures/BB_p_Forest.tiff');
cc = imread('../Figures/CC_p_Forest.tiff');
dd = imread('../Figures/DD_p_Forest.tiff');
isis = imread('../Figures/ISIS_p_Forest.tiff');
oo = imread('../Figures/OO_p_Forest.tiff');
all = imread('../Figures/4x4_Forest.tiff');

inds = 2500:3350;
gap = ones(5354, 200, 3)*255;

combined = [all(:, 1:1780, :) gap bb(:, inds, :) gap aa(:, inds, :) gap cc(:, inds, :) ...
    gap dd(:, inds, :) gap isis(:, inds, :) gap oo(:, inds, :) gap gap];

REModelLocations = [ ...
    1827:1973 ...
    2542:2688 ...
    3183:3329 ...
    3969:4115 ...
    4184:4330 ...
    4685:4831];
OverallRELocation = 4826:4972;

combined(REModelLocations,:,:) = combined(REModelLocations,:,:)*0.75;
combined(OverallRELocation,:,:) = combined(OverallRELocation,:,:)*0.55;

combined = combined(308:end, 65:8370, :);
combined(1:200, 1:5000, :) = 255;

image(combined)

imwrite(combined, '../Figures/Combined_p_underlay.tiff')

szC = size(combined);
dpi = 600;
sz = [0 0 szC([2 1])];

f = figure(...
  'PaperUnits','inches',...
  'PaperPosition', sz/dpi,...
  'PaperPositionMode','manual',...
  'Visible',  'off');
axes('Position', [0 0 1 1])
image(255*ones(size(combined)))

[~, summaryEffects] = xlsread('../Figures/ContingencyTables.xlsm', 'Data for forests');
for i = 1:size(summaryEffects, 2)
    text_p = summaryEffects{1, i};
    num_p = round(str2double(text_p), 3, 'significant');
    summaryEffects{1, i} = num2str(num_p);
    
    text_ci = summaryEffects{2, i};
    if contains(text_ci, 'NaN')
        summaryEffects{2, i} = '';
    else
        text_ci = strsplit(text_ci, ', ');
        num_cil = round(str2double(text_ci{1}(2:end)), 3, 'significant');
        num_ciu = round(str2double(text_ci{2}(1:end-1)), 3, 'significant');
        summaryEffects{2, i} = ['[' num2str(num_cil) ', ' num2str(num_ciu) ']'];
    end
end

SummaryLocationsX = [2670 3726 4775 5823 6876 7929]'*ones(1,7);
SummaryLocationsX = SummaryLocationsX(:);
SummaryLocationsY = [[REModelLocations(1:147:end)-275 OverallRELocation(1)-275]'*ones(1,6)]';
SummaryLocationsY = SummaryLocationsY(:);

text(SummaryLocationsX(1:36), SummaryLocationsY(1:36), summaryEffects(1,1:36), ...
    'FontSize', fontSize)
text(SummaryLocationsX(37:end), SummaryLocationsY(37:end), summaryEffects(1,37:end), ...
    'FontSize', fontSize, 'FontWeight', 'bold')

text(SummaryLocationsX(1:36), SummaryLocationsY(1:36)+70, summaryEffects(2,1:36), ...
    'FontSize', fontSize)
text(SummaryLocationsX(37:end), SummaryLocationsY(37:end)+70, summaryEffects(2,37:end), ...
    'FontSize', fontSize, 'FontWeight', 'bold')

text(1620, 100, 'N', 'FontSize', 9, 'FontWeight', 'bold')
text(2220, 100, 'Secure', 'FontSize', 9, 'FontWeight', 'bold')
text(3245, 100, 'Avoidant', 'FontSize', 9, 'FontWeight', 'bold')
text(4255, 100, 'Ambivalent', 'FontSize', 9, 'FontWeight', 'bold')
text(5270, 100, 'Disorganized', 'FontSize', 9, 'FontWeight', 'bold')
text(6400, 100, 'Insecure', 'FontSize', 9, 'FontWeight', 'bold')
text(7415, 100, 'Organized', 'FontSize', 9, 'FontWeight', 'bold')

% axis tight
axis off
% axis equal
print(f, '-r600','-dtiff','../Figures/Combined_p_overlay.tiff');