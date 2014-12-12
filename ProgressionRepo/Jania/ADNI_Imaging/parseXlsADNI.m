function parseXlsADNI

%import the DNI Clinical data
[num,txt,raw] = xlsread('/Users/jaghajan/Documents/Experiments/MAPER_DATA/adni_pdxconv_2011-06-02.xls');

%read the second column of num which is the patent ID i.e. RID
patientID=num(:,2);
%read the fifth column of num wihch is the current diagnosis i.e. DXCURRENT
%if DXCURREN == 1 then healthy control patient
dxcurren=num(:,5);

%read the thr=ird column which is the visitID to pick out only the baseline
%Note exclude the first row as it is the header
visitID=txt(2:end,3);
visitBL=find(ismember(visitID, 'bl')==1);

%store the diagnosis info for baseline
patientIDBaseline=patientID(visitBL);
dxcurrenBaseline=dxcurren(visitBL);
save('/Users/jaghajan/Documents/Experiments/MAPER_DATA/DiagnosisInfo.mat','patientIDBaseline','dxcurrenBaseline');

%select the control patients
indx=find(dxcurrenBaseline==1);
controlPatientID=patientIDBaseline(indx);

save('/Users/jaghajan/Documents/Experiments/MAPER_DATA/ControlPatientID.mat','controlPatientID');
