function [all_pop_matrix] = calc_variance_diagrams(samples,max_seq_grad)

NR_BIOMK = size(max_seq_grad,2);
all_pop_matrix = zeros(NR_BIOMK, NR_BIOMK); % all_pop_matrix(event, position);

for s=1:size(samples,1)
  for biomk = 1:NR_BIOMK
    all_pop_matrix(biomk, samples(s,biomk)) = all_pop_matrix(biomk, samples(s,biomk)) + 1;
  end
end

%permute the matrix using the characteristic/max_likelihood ordering
 
all_pop_matrix = all_pop_matrix(max_seq_grad,:);

end