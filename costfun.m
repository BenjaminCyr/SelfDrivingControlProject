function [J, dJ] = costfun(z, nsteps, final_state, weights)
    
    x_goal = final_state(1);
    y_goal = final_state(3);
    
    x = z(1:6:(nsteps)*6);
    y = z(3:6:(nsteps)*6);
    Fx = z((nsteps)*6+1:2:(nsteps)*8-2);
    df = z((nsteps)*6+2:2:(nsteps)*8-2);
    T = z(end);

    cost_values = [(x(end) - x_goal)^2 + (y(end) - y_goal)^2; sum(Fx.^2); sum(df.^2); T^2];
    J = weights*cost_values;
    
    
    dJ = zeros(length(z),1);
    
    dJ(nsteps*6-5) = weights(1)*2*(x(end) - x_goal);
    dJ(nsteps*6-3) = weights(1)*2*(y(end) - y_goal);
    dJ((nsteps)*6+1:2:(nsteps)*8-2) = weights(2)*2*Fx;
    dJ((nsteps)*6+2:2:(nsteps)*8-2) = weights(3)*2*df;
    dJ(end) = 2*weights(4)*T;
end
