function [J, dJ] = costfun(z, nsteps, final_state, weights)
%     [a, b, c] = get_line(clinePt, final_state(1:2:3));
    
    x_goal = final_state(1);
    y_goal = final_state(3);
    
    x = z(1:6:(nsteps)*6);
    y = z(3:6:(nsteps)*6);
    Fx = z((nsteps)*6+1:2:(nsteps)*8-2);
    df = z((nsteps)*6+2:2:(nsteps)*8-2);
    
%     dist = (a*x(1:end-1) + b*y(1:end-1) + c)/sqrt(a^2 + b^2);
%     sum(dist.^2);

    J = sum(weights*[(x(end) - x_goal)^2 + (y(end) - y_goal)^2; sum(Fx.^2); sum(df.^2)]);
    
    
    dJ = zeros(8*nsteps-2,1);
%     dJ(1:6:(nsteps-1)*6) = weights(1)*(2*a^2*x(1:end-1) + 2*a*b*y(1:end-1) + 2*a*c)/(a^2 + b^2);
%     dJ(3:6:(nsteps-1)*6) = weights(1)*(2*b^2*y(1:end-1) + 2*a*b*x(1:end-1) + 2*b*c)/(a^2 + b^2);
    
    dJ(nsteps*6-5) = weights(1)*2*(x(end) - x_goal);
    dJ(nsteps*6-3) = weights(1)*2*(y(end) - y_goal);
    dJ((nsteps)*6+1:2:(nsteps)*8-2) = weights(2)*2*Fx;
    dJ((nsteps)*6+2:2:(nsteps)*8-2) = weights(3)*2*df;
end
