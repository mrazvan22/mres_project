div_size_list=[0.1,0.3,0.5,0.75,0.8,0.85,0.9,0.95];
version_likelihood=8; 
n_repeat=10;
for i=1:length(div_size_list)
div_size=div_size_list(i)
robustnessTest2Demo(div_size,version_likelihood,n_repeat)
pause
end