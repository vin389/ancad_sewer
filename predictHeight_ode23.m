function ht_t56_predict = ...
    predictHeight_ode23(t5, t6, rain_mm, ht_m, lambda, hb, gamma)
    % trim the whole data (rain_mm and ht_m) to trimmed data (i.e., rain_t56
    % and ht_t56)
    t56 = double((t5:t6) - t5); t56 = t56(:);
    ht_t56 = ht_m(t5:t6); ht_t56 = ht_t56(:);   
    rain_t56 = rain_mm(t5:t6); rain_t56 = rain_t56(:);
    % debug for development 
    %     calculate predicted ht_t56 (ht_t56_predict) based on lambda, hb, gamma
    theDebug = 0;
    if (theDebug >= 1)
        ht_t56_predict_forward = zeros(size(ht_t56)); % the _cd stands for forward explicit method
        ht_t56_predict_forward(1) = ht_t56(1);
        for i = 2: size(ht_t56_predict_forward(:), 1)
            % predict water level step by step
            ht_t56_predict_forward(i) = ht_t56_predict_forward(i - 1) + ...
                gamma * rain_t56(i - 1) - lambda * sqrt(ht_t56_predict_forward(i - 1) - hb);
            % minimum of water level is zero
            if (ht_t56_predict_forward(i) <= max(0., hb))
                ht_t56_predict_forward(i) = max(0., hb); 
            end
        end
    end
    % predict ht_t56 using ode23 
    ht_t56_predict_ode23 = zeros(size(ht_t56)); 
    % ode23 usage: [t,y] = ode23(odefun,tspan,y0)
    tspan = [t56(1) t56(end)];
    [t_ode23, ht_ode23] = ode23(@(t, ht_t) predictHeight_ode23_odefun( ...
        t, ht_t, t56, rain_t56, lambda, hb, gamma), tspan, ht_t56(1));
    ht_t56_predict_ode23 = interp1(t_ode23, ht_ode23, t56);
    % debug for development
    %     compare ht_t56_predict_forward and ht_ode23
    if (theDebug >= 1) 
        figure; plot(t56, ht_t56_predict_forward, ...
                     t_ode23, ht_ode23, ...
                     t56, ht_t56_predict_ode23); grid('on'); 
        legend('Finite diff. forward', 'ODE 23', 'ODE 23 at integer t');
    end
    % return
    ht_t56_predict = ht_t56_predict_ode23;
end
