function [ht_t56_predict] = ...
   predictHeight(t5, t6, rain_mm, ht_m, lambda, hb, gamma)
   % trim the whole data (rain_mm and ht_m) to trimmed data (i.e., rain_t56
   % and ht_t56)
   t56 = double((t5:t6) - t5); t56 = t56(:);
   ht_t56 = ht_m(t5:t6); ht_t56 = ht_t56(:);   
   rain_t56 = rain_mm(t5:t6); rain_t56 = rain_t56(:);
   % calculate predicted ht_t56 (ht_t56_predict) based on lambda, hb, gamma
   ht_t56_predict = zeros(size(ht_t56));
   ht_t56_predict(1) = ht_t56(1);
   for i = 2: size(ht_t56_predict(:), 1)
       % predict water level step by step
       ht_t56_predict(i) = ht_t56_predict(i - 1) + ...
           gamma * rain_t56(i - 1) - lambda * sqrt(ht_t56_predict(i - 1) - hb);
       % minimum of water level is zero
       if (ht_t56_predict(i) <= max(0., hb))
           ht_t56_predict(i) = max(0., hb); 
       end
   end
end
