function [percentile_vec] = generate_percentiles(toanalyze_vec)
%Returns the percentile vector from a numeric vector
[~, sort_i] = sort(toanalyze_vec);
[~, vec_rank] = sort(sort_i);

percentile_vec = vec_rank/length(vec_rank);


end

