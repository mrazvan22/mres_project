function h = plotVarMatrix(mat, ordering)

labels =  {'T-tau', 'Abeta', 'P-tau', 'Ventricles','Hippocampus',...
  'WholeBrain','Entorhinal','Fusiform','Mid Temporal','ADAS-Cog','MMSE',...
  'RAVLT','Brain Atrophy','Hippo. atrophy'};
h = figure
imagesc(1-mat);
colormap(gray);
ax = gca;
ax.XTick = 1:14;
ax.YTick = 1:14;
ax.YTickLabel = labels(ordering);
end