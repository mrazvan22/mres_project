function [max_seq,max_lik,final_sequences,final_lik] = char_seq_grad_asc(X, mu_mix, sigma_mix, pi_mix)

[NR_SUBJECTS,NR_BIOMK] = size(X);


[nr_biomk2,~] = size(mu_mix);
[nr_biomk3,~] = size(sigma_mix);

assert(NR_BIOMK == nr_biomk2 && NR_BIOMK == nr_biomk3);


NR_ITERATIONS = 2000;
ROUNDS = 10; % start from 10 different initialisation points in order to ensure we find the global maximum

s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

final_sequences = zeros(ROUNDS, NR_BIOMK);
final_lik = zeros(ROUNDS, 1);
last_changed = zeros(ROUNDS,1);

for r=1:ROUNDS
  r
  curr_seq = randperm(NR_BIOMK);
  old_likelihood = calc_likelihood(X, curr_seq, mu_mix, sigma_mix, pi_mix);
  for i=1:NR_ITERATIONS
      p1 = ceil(rand * (NR_BIOMK - 1));
      p2 = p1+1;
      %p2 = ceil(rand * NR_BIOMK);
      %while (p1 == p2)
      %    p2 = ceil(rand * NR_BIOMK);
      %end
      new_seq = curr_seq;

      % swap events
      tmp = new_seq(p1);
      new_seq(p1) = new_seq(p2);
      new_seq(p2) = tmp;

      new_likelihood = calc_likelihood(X, new_seq, mu_mix, sigma_mix, pi_mix);
      if(new_likelihood >= old_likelihood)
          curr_seq = new_seq;
          old_likelihood = new_likelihood;
          last_changed(r) = i;
      end

%       curr_seq

  end
  final_sequences(r,:) = curr_seq;
  final_lik(r) = old_likelihood;
end

[max_lik, maxIndex] = max(final_lik);
max_seq = final_sequences(maxIndex,:);

end