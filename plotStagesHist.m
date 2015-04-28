function h = plotStagesHist(patientStages, EBMdxBL)

h = figure

controlIndices = find(EBMdxBL == 1);
mciIndices = find(EBMdxBL == 2);
adIndices = find(EBMdxBL == 3);

controlStages = patientStages(controlIndices);
mciStages = patientStages(mciIndices);
adStages = patientStages(adIndices);

nrBins = 15;
contHist = hist(controlStages,nrBins);
mciHist = hist(mciStages,nrBins);
adHist = hist(adStages,nrBins);

% normalise the values
contHist = contHist / sum(contHist);
mciHist = mciHist / sum(mciHist);
adHist = adHist / sum(adHist);


b = bar([contHist; mciHist; adHist]');
b(1).FaceColor = 'blue';
b(2).FaceColor = 'k';
b(3).FaceColor = 'red';


legend('control', 'MCI', 'AD');
end