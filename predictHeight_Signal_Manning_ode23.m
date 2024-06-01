function ht_t910_predict = ...
    predictHeight_Signal_Manning_ode23(t9, t10, rain_mm, ht_m, lambda, ...
    hb, gamma, diameter, knss)
    % trim the whole data (rain_mm and ht_m) to trimmed data (i.e.,
    % rain_t910 and ht_t910)
    t910 = double((t9:t10) - t9); t910 = t910(:);
    ht_t910 = ht_m(t9:t10); ht_t910 = ht_t910(:);
    rain_t910 = rain_mm(t9:t10); rain_t910 = rain_t910(:);
    % predict ht_t56 using ode23 
    ht_t910_predict_ode23 = zeros(size(ht_t910)); 
    % ode23 usage: [t,y] = ode23(odefun,tspan,y0)
    tspan = [t910(1) t910(end)];
    [t_ode23, ht_ode23] = ode23(@(t, ht_t) ...
        predictHeight_Signal_Manning_ode23_odefun( ...
        t, ht_t, t910, rain_t910, lambda, hb, gamma, diameter, knss), ...
        tspan, ht_t910(1));
    ht_t910_predict_ode23 = interp1(t_ode23, ht_ode23, t910);
    % debug for development
    %     compare ht_t56_predict_forward and ht_ode23
    theDebug = 0; 
    if (theDebug >= 1) 
        figure; plot(t56, ht_t56_predict_forward, ...
                     t_ode23, ht_ode23, ...
                     t56, ht_t56_predict_ode23); grid('on'); 
        legend('Finite diff. forward', 'ODE 23', 'ODE 23 at integer t');
    end
    % return
    ht_t910_predict = ht_t910_predict_ode23;
end
