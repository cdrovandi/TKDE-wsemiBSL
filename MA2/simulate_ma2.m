function y = simulate_ma2(theta,T,reps)
    % function to simulate MA2 model (can simulate reps times uisng vectorised code)
    e = randn(T+2,reps);
    y = e(3:end,:) + theta(1)*e(2:end-1,:) + theta(2)*e(1:end-2,:);
end


