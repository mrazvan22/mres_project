function char_sequence = calc_char_sequence(samples)

biomk_avg_ordering = mean(samples);
[~, char_sequence] = sort(biomk_avg_ordering);