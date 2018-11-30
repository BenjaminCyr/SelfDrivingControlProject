function [g,h,dg,dh]=nonlcon(z, nsteps, initial_state, T)
    dt = T/(nsteps - 1);
    
    %load('TestTrack.mat');
    %left = TestTrack.br;
    %right = TestTrack.bl;
    
    
    x = z(1:6:(nsteps)*6);
    u = z(2:6:(nsteps)*6);
    y = z(3:6:(nsteps)*6);
    v = z(4:6:(nsteps)*6);
    psi = z(5:6:(nsteps)*6);
    r = z(6:6:(nsteps)*6);
    
    Fx = z((nsteps)*6+1:2:(nsteps)*8-2);
    delta = z((nsteps)*6+2:2:(nsteps)*8-2);
    
    g = zeros(nsteps, 1); 

    dg = zeros(length(z), length(g));
   
    % size of h must be ((no. of time steps * no. of states) x 1)
    h = zeros(nsteps*6, 1);
    h(1:6) = z(1:6) - initial_state';
    
    % size of dh must be Transpose((no. of time steps * no. of states) x no. of elements in 'z') ;
    dh = zeros(length(z), length(h));
    dh(1:6,1:6) = diag(ones(1,6));
    
    for i = 2:nsteps-1
        base_i = 6*i-5;
       dx_i_1 = bike_odefun([x(i-1) u(i-1) y(i-1) v(i-1) psi(i-1) r(i-1)], [Fx(i-1) delta(i-1)]);
       dx_i = bike_odefun([x(i) u(i) y(i) v(i) psi(i) r(i)], [Fx(i) delta(i)]);

       x_i = [x(i) u(i) y(i) v(i) psi(i) r(i)]';
       x_i_1 = [x(i-1) u(i-1) y(i-1) v(i-1) psi(i-1) r(i-1)]';
       
       h(base_i:base_i+5) = x_i - x_i_1 - 0.5*dt*(dx_i + dx_i_1);
%        h(base_i:base_i+5) = x_i - x_i_1 - dt*dx_i_1;
        
       %x_i
       dh(base_i:base_i+5,base_i:base_i+5) = diag(ones(1,6));
       %x_i-1
       dh(base_i-6:base_i-1, base_i:base_i+5) = diag(-ones(1,6));
                                
       %x_dot_i-1
       [A_i_1, B_i_1] = calculate_jacobian([x(i-1) u(i-1) y(i-1) v(i-1) psi(i-1) r(i-1)], [Fx(i-1) delta(i-1)]);
       [A_i, B_i] = calculate_jacobian([x(i) u(i) y(i) v(i) psi(i) r(i)], [Fx(i) delta(i)]);
       
       dh(base_i-6:base_i-1, base_i:base_i+5) = dh(base_i-6:base_i-1, base_i:base_i+5) - 0.5*dt*(A_i_1');
       dh(2*i+(nsteps)*6-3:2*i+(nsteps)*6-2, base_i:base_i+5) = -0.5*dt*(B_i_1');
       
       dh(base_i:base_i+5, base_i:base_i+5) = dh(base_i:base_i+5, base_i:base_i+5) - 0.5*dt*(A_i');
       dh(2*i+(nsteps)*6-1:2*i+(nsteps)*6, base_i:base_i+5) = -0.5*dt*(B_i'); 
       
%        dh(end, base_i:base_i+5) = -0.5*(dx_i + dx_i_1)./(nsteps-1);
%        dh(end, base_i:base_i+5) = -dx_i_1./(nsteps-1);
    end
    
     i = nsteps;
     base_i = 6*i-5;
   dx_i_1 = bike_odefun([x(i-1) u(i-1) y(i-1) v(i-1) psi(i-1) r(i-1)], [Fx(i-1) delta(i-1)]);
   x_i = [x(i) u(i) y(i) v(i) psi(i) r(i)]';
   x_i_1 = [x(i-1) u(i-1) y(i-1) v(i-1) psi(i-1) r(i-1)]';
   h(base_i:base_i+5) = x_i - x_i_1 - dt*dx_i_1;

   %x_i
   dh(base_i:base_i+5,base_i:base_i+5) = diag(ones(1,6));
   %x_i-1
   dh(base_i-6:base_i-1, base_i:base_i+5) = diag(-ones(1,6)); 

   %x_dot_i-1
   [A_i_1, B_i_1] = calculate_jacobian([x(i-1) u(i-1) y(i-1) v(i-1) psi(i-1) r(i-1)], [Fx(i-1) delta(i-1)]);

   dh(base_i-6:base_i-1, base_i:base_i+5) = dh(base_i-6:base_i-1, base_i:base_i+5) - dt*(A_i_1');
   dh(2*i+(nsteps)*6-3:2*i+(nsteps)*6-2, base_i:base_i+5) = -dt*(B_i_1');
%    dh(end, base_i:base_i+5) = -dx_i_1./(nsteps-1);
end