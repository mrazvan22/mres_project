function [mu_mix, sigma_mix, pi_mix] = calc_gaussian_parameters(EBMdataBL, EBMdxBL)


[nr_patients, nr_biomarkers] = size(EBMdataBL);

control_indices = find(EBMdxBL == 1);
patient_indices = find(EBMdxBL > 1);

mu_mix = zeros(nr_biomarkers, 2);
sigma_mix = zeros(nr_biomarkers, 2);
pi_mix = zeros(nr_biomarkers, 2);

for biomk=1:nr_biomarkers
   control_levels = EBMdataBL(control_indices, biomk);
   patient_levels = EBMdataBL(patient_indices, biomk);
    
   mu_control = mean(control_levels);
   sigma_control = std(control_levels);
   
   mu_patient = mean(patient_levels);
   sigma_patient = std(patient_levels);
   
   %pXgEnE(:, biomk, 2) = normpdf(control_levels, mu_control, sigma_control);
   %pXgEnE(:, biomk, 1) = normpdf(patient_levels, mu_patient, sigma_patient);
   
   
   all_levels = EBMdataBL(:, biomk);
   % fit two gaussians on all the data

   [mu_mix(biomk,:), sigma_mix(biomk,:), pi_mix(biomk,:)]  = ...
       em_mix(all_levels, [mu_control, mu_patient], [sigma_control, sigma_patient])
   

   minX = min(all_levels);
   maxX = max(all_levels);
   X = minX:(maxX - minX)/200:maxX;
   Y = eval_mix_gaussians(X, mu_mix(biomk,:), sigma_mix(biomk,:), pi_mix(biomk,:));
   %Y = eval_mix_gaussians(X, [mu_control, mu_patient], [sigma_control, sigma_patient], pi_mix(biomk,:));
   
   
   clf;
   [f1,x1] = hist(patient_levels);
   bar(x1,f1, 'FaceColor',[0.8,0.2,0.2]);
   %h = findobj(gca,'Type','patch');
   %h.FaceColor = [0.8 0.2 0.2];
   %h.EdgeColor = 'b';
   hold on
   [f2,x2] = hist(control_levels);
   bar(x2,f2,'b');
   maxY = get(gca,'ylim');
   Y = Y * (maxY(2) * 0.66/max(Y));
   %hist(control_levels);
   hold on

   plot(X, Y,'y');

   %hist(control_levels,50)
   
   if(biomk == 14)
        display('asfda')
   end
   
end

   save('sporadic_params.mat', 'mu_mix', 'sigma_mix', 'pi_mix'); 
   

end