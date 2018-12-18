function U = ROB599_ControlsProject_part2_Team35(TestTrack, Xobs) 
    timeout = tic;
    centroids = reshape(mean(cell2mat(Xobs), 1), 2, length(Xobs))';
    centerline = TestTrack.cline;
    leftline = TestTrack.bl;
    rightline = TestTrack.br;
    
    final_point_left = [leftline(:,end)+[20;50]];
    final_point_right = [rightline(:,end)+[20;50]];
    final_point = [mean([final_point_left(1) final_point_right(1)]);
                   mean([final_point_left(2) final_point_right(2)])];
    centerline = [centerline final_point];
    leftline = [leftline final_point];
    rightline = [rightline final_point];

    w = [0.6 0.4];
    x_l = w*[centerline(1,:); leftline(1,:)];
    y_l = w*[centerline(2,:); leftline(2,:)];
    leftline = [x_l; y_l];
    
    x_r = w*[centerline(1,:); rightline(1,:)];
    y_r = w*[centerline(2,:); rightline(2,:)];
    rightline = [x_r; y_r];
    
    lines = [leftline; centerline; rightline];

    traj_total = generate_trajectory(lines, centroids);
%     load('test.mat', 'traj_total');
    
    U = MPC_control(traj_total);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% MPC Control
    function U_out = MPC_control(traj_total)
        ub_df = 0.5;
        lb_df = -0.5;

        ub_Fx = 2500;
        lb_Fx = -5000;

        dt=0.1;
        dt_final = 0.01;
        num_steps = dt/dt_final;

        options = optimoptions('quadprog', 'Display', 'none');

        T = traj_total(:,1);
        %time span
        tspan = [0:dt:T(end)]';

        start_point = 1;
        end_point = length(tspan);

        %Interpolate;
        Z_ref=interp1(T, traj_total(:,2:9),tspan);
        Y_ref = Z_ref(:,1:6);
        U_ref = Z_ref(:, 7:8);

        state_c = @(t) Z_ref(round(t/dt+1), 1:6);
        input_c = @(t) Z_ref(round(t/dt+1), 7:8);

        A_c = @(t) calculate_A(state_c(t), input_c(t));
        B_c = @(t) calculate_B(state_c(t), input_c(t));

        A=@(i) eye(6) + dt*A_c((i-1)*dt);
        B=@(i) dt*B_c((i-1)*dt);

        %% 2.2 Number of decision variables for colocation method
        npred=10;

        %% 2.5 simulate controller working from initial condition [0.25;-0.25;-0.1]
        %use ode45 to between inputs
        % we begin by defining the cost function
        Q = diag([10 1 10 1 1 1]);
        R = diag([.000001 1]);

        H = zeros(6*(npred+1) + 2*npred);
        c = zeros(6*(npred+1) + 2*npred,1);

        H((6*npred+1):(6*npred+6),(6*(npred)+1):(6*npred+6)) = Q;
        for i = 1:npred
            H((6*(i-1)+1):(6*(i-1)+6),(6*(i-1)+1):(6*(i-1)+6)) = Q;
            H((2*(i-1)+6*(npred+1)+1):(2*i+6*(npred+1)), (2*(i-1)+6*(npred+1)+1):(2*i+6*(npred+1))) = R;
        end

        Y = zeros(length(tspan), 6);
        U_out = zeros(length(tspan), 2);
        Y(start_point,:) = Y_ref(start_point,:);

        last_U = [U_ref(start_point,2) U_ref(start_point,1)];

        for i = start_point:end_point %length(tspan)-1 
           error = Y(i, :)-Y_ref(i,:);
           if abs(error(1)) > 15 || toc(timeout) > 14.5*60
              break; 
           end
%            fprintf('Error %d:  [%f %f %f %f %f %f]\n', i, error);
           window = min(npred, length(tspan)-i);
           if window > 1
               [Aeq, beq] = eq_cons(i,A,B,window, error);
               [Lb, Ub] = bound_cons(i, U_ref, [lb_Fx ub_Fx; lb_df ub_df], window);

               H_section = H((npred-window)*6+1 : 6*npred+6+2*window, (npred-window)*6+1 : 6*npred+6+2*window);
               optimal_z = quadprog( H_section, c(1:8*window+6), [], [], Aeq, beq, Lb, Ub, [], options);
           end

           y_mpc = optimal_z(1:6*(window+1));
           u_mpc = optimal_z(6*(window+1)+1:6*(window+1)+2);
           current_U = [U_ref(i,2)+u_mpc(2) U_ref(i, 1)+u_mpc(1)];
           u = interp1([0 dt], [last_U; current_U], 0:dt_final:dt);
           %use ode45 to simulate nonlinear system, f, forward 1 timestep
           ytemp = forwardIntegrateControlInput2(u, Y(i,:));
        %    [~, ytemp] = ode45(@(t, x) bike_odefun(x, current_U), [0 dt], Y(i,:));
           Y(i+1,:) = ytemp(end,:);
        %    
        %    plot(Y(i:i+1, 1),Y(i:i+1, 3), 'r')
        %    hold on

           U_out(i*num_steps-num_steps+1:i*num_steps,:) = u(1:end-1, :);
           last_U = current_U;
        end
    end

    function [Aeq,beq]=eq_cons(idx,A,B,npred,x0)
        %build matrix for A_i*x_i+B_i*u_i-x_{i+1}=0
        %in the form Aeq*z=beq
        %initial_idx specifies the time index of initial condition from the reference trajectory 
        %A and B are function handles above

        n = size(A(1), 2); %size of state vector
        m = size(B(1), 2); % size of input vector

        %Aeq*[x_0, y_0, psi_0, ... ,u_0, delta_0, ...] =beq
        Aeq = zeros(n*(npred+1), n*(npred+1)+m*npred);
        beq = zeros(n*(npred+1), 1);

        Aeq( 1:n, 1:n ) = eye( n ); 
        beq( 1:n ) = x0';
        for i = 1:npred-1
            Aeq((i*n+1):((i+1)*n), (i*n+1):((i+1)*n)) = eye(n);
            Aeq((i*n+1):((i+1)*n), ((i-1)*n+1):(i*n)) = -A(idx+i-1);
            Aeq((i*n+1):((i+1)*n), ((i-1)*m+(npred+1)*n+1):(i*m+(npred+1)*n)) = -B(idx+i-1);
        end
    end

    function [Lb,Ub]=bound_cons(idx, U_ref, input_range, npred)
        %initial_idx is the index along uref the initial condition is at
        Lb = -Inf(6*(npred+1)+2*npred, 1);
        Ub = Inf(6*(npred+1)+2*npred, 1);

        for i = 1:npred
            Lb(2*i+6*(npred+1)-1:2*i+6*(npred+1)) = input_range(:,1) - U_ref(idx+i-1,:)';
            Ub(2*i+6*(npred+1)-1:2*i+6*(npred+1)) = input_range(:,2) - U_ref(idx+i-1,:)';
        end
    end
    
    
    %FMINCON trajectory generation
    function traj_total = generate_trajectory(lines, centroids)
        cline = lines(3:4,:);
        
        radius = 4;
        
        nsteps = 15+1;
        min_T = 0.01*nsteps;
        max_T = 10;

        highest_x = 1500;
        lowest_x = 0;

        highest_y = 850;
        lowest_y = -200;

        ub_delta = 0.5;
        lb_delta = -0.5;

        ub_Fx = 2500;
        lb_Fx = -5000;

        % weights = [1 0.1 1 0 0 0 0 0 nsteps];
        weights = [1 .00001 1 5];


        ub = [repmat([highest_x Inf highest_y Inf Inf Inf], 1, nsteps) repmat([ub_Fx ub_delta], 1, nsteps-1) max_T];

        lb = [repmat([lowest_x -Inf lowest_y -Inf -Inf -Inf], 1, nsteps) repmat([lb_Fx lb_delta], 1, nsteps-1) min_T];

        % %1.4
        % %%%%%%%%%%%%%%% no need to change these lines  %%%%%%%%%%%%%%%%%%%%%%
        options = optimoptions('fmincon', 'ConstraintTolerance', 1e-6, 'SpecifyConstraintGradient',true,...
                               'SpecifyObjectiveGradient',true,...
                               'MaxFunctionEvaluations', 3000, 'MaxIterations', 10000, 'Display','none',...
                               'CheckGradients',false);
        %                          'Display','iter')

        is = [287 5 -176 0 2 0];

        fs = [0 0 0 0 0 0];

        traj_total  = [];
        last_time = 0;
        index = 1;
        i_iter = 1;
        while is(3) < cline(2,end-1)
            [close_point, index] = closest_point([is(1) is(3)], lines, index, centroids);
            fs = [close_point(1) 0 close_point(2) 0 0 0];
            
%             plot(is(1),is(3), 'ro')
%             plot(fs(1),fs(3), 'go')

            x0 = linspace(is(1), fs(1), nsteps);
            y0 = linspace(is(3), fs(3), nsteps);

            z0 = zeros(nsteps*8-2+1, 1);
            z0(1:nsteps*6) = repmat(is, 1, nsteps);
            z0(1:6:nsteps*6) = x0;
            z0(3:6:nsteps*6) = y0;

            z0(end) = 2;

            nc=@(z) nonlcon(z, nsteps, is, centroids, radius);
            cf=@(z) costfun(z, nsteps, fs, weights);

%             fprintf('%d:\t', i_iter);
            [z,fval,~,output] = fmincon(cf,z0,[],[],[],[],lb',ub',nc,options);
%             fprintf('fval = %e\tFeasibility = %e\n', fval, output.constrviolation);
            if output.constrviolation > 1e-1 || toc(timeout) > 10*60
               break; 
            end

            Y0=reshape(z(1:6*nsteps),6,nsteps)';
            U0=reshape(z(6*nsteps+1:8*nsteps-2),2,nsteps-1);

            T_iter = z(end);
            dt_iter = T_iter/(nsteps - 1);


            dt_total = linspace(last_time, last_time+T_iter, nsteps);
            traj_total = [traj_total;dt_total(1:end-1)' Y0(1:end-1, :) U0']; 
            last_time = last_time+T_iter;
            
            is = Y0(end,:);
            i_iter = i_iter + 1;
            
%             plot(Y0(:,1),Y0(:,3), 'c')
%             hold on
        end

%         traj_total(end,1)
    end

    function [pt, current_index] = closest_point(pt_in, lines, current_index, centroids)
        lline = lines(1:2,:);
        cline = lines(3:4,:);
        rline = lines(5:6,:);
        
        distance = sqrt((pt_in(1) - cline(1,current_index))^2 + (pt_in(2) - cline(2,current_index))^2);
        while distance < 20 && current_index < length(cline) %|| min(dist_centroids) < radius   
           current_index = current_index+1;
           distance = sqrt((pt_in(1) - cline(1,current_index))^2 + (pt_in(2) - cline(2,current_index))^2);
        end
        pt = cline(:,current_index);
        
        dist_centroids = sqrt((centroids(:,1) - pt(1)).^2 + (centroids(:,2) - pt(2)).^2);
        if min(dist_centroids) < 40
            pt_l = lline(:,current_index);
            dist_centroids_left = sqrt((centroids(:,1) - pt_l(1)).^2 + (centroids(:,2) - pt_l(2)).^2);
            pt_r = rline(:,current_index);
            dist_centroids_right = sqrt((centroids(:,1) - pt_r(1)).^2 + (centroids(:,2) - pt_r(2)).^2);
            if min(dist_centroids_left) > min(dist_centroids_right)
                pt = pt_l;
            else
                pt = pt_r;
            end
        end
    end
    
    function [g,h,dg,dh]=nonlcon(z, nsteps, initial_state, centroids, radius)
        T = z(end);
        dt = T/(nsteps - 1);

        x = z(1:6:(nsteps)*6);
        u = z(2:6:(nsteps)*6);
        y = z(3:6:(nsteps)*6);
        v = z(4:6:(nsteps)*6);
        psi = z(5:6:(nsteps)*6);
        r = z(6:6:(nsteps)*6);

        Fx = z((nsteps)*6+1:2:(nsteps)*8-2);
        delta = z((nsteps)*6+2:2:(nsteps)*8-2);

        %(nsteps*num_centroids) x 1
        g = zeros(nsteps*size(centroids, 1), 1);
        for i_obj = 0:size(centroids, 1)-1
            g(nsteps*i_obj+1:nsteps*(i_obj+1)) = radius^2 - (x - centroids(i_obj+1,1)).^2 - (y - centroids(i_obj+1,2)).^2; 
        end

        %z x (nsteps*num_centroids)
        dg = zeros(length(z), length(g));
        for i_obj = 0:size(centroids, 1)-1
            dg(1:6:6*nsteps, i_obj*nsteps+1:nsteps*(i_obj+1)) = diag(-2*(x - centroids(i_obj+1,1)));
            dg(3:6:6*nsteps, i_obj*nsteps+1:nsteps*(i_obj+1)) = diag(-2*(y - centroids(i_obj+1,2)));
        end

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

           %T
           dh(end, base_i:base_i+5) = -0.5*(dx_i + dx_i_1)./(nsteps-1);
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
       dh(end, base_i:base_i+5) = -dx_i_1./(nsteps-1);
    end
    
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

    function [A, B] = calculate_jacobian(state, input)
        A = calculate_A(state, input);
        B = calculate_B(state, input);
    end

    function A = calculate_A(state, input)
        Nw=2;
        f=0.01;
        Iz=2667;
        a=1.35;
        b=1.45;
        By=0.27;
        Cy=1.2;
        Dy=0.7;
        Ey=-1.6;
        Shy=0;
        Svy=0;
        m=1400;
        g=9.806;

        X = state(1); u = state(2);
        Y = state(3); v = state(4);
        psi = state(5); r = state(6);

        if u == 0
           u = 0.001; 
        end
        %generate input functions
        Fx=input(1);
        delta=input(2);

        %slip angle functions in degrees
        a_f=rad2deg(delta-atan2(v+a*r,u));
        a_r=rad2deg(-atan2((v-b*r),u));

        %Alpha derivatives
        da_fdu = rad2deg((v+a*r)/((v+a*r)^2 + u^2));
        da_fdv = rad2deg(-u/((v+a*r)^2 + u^2));
        da_fdr = rad2deg(-a*u/((v+a*r)^2 + u^2));

        da_rdu = rad2deg((v-b*r)/((v-b*r)^2 + u^2));
        da_rdv = rad2deg(-u/((v-b*r)^2 + u^2));
        da_rdr = rad2deg(b*u/((v-b*r)^2 + u^2));

        %Nonlinear Tire Dynamics
        phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
        phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

        %Phi derivatives
        dphi_yfda_f = (1-Ey) + Ey/(By^2*(a_f+Shy)^2 + 1);
        dphi_yrda_r = (1-Ey) + Ey/(By^2*(a_r+Shy)^2 + 1);


        %Lateral Forces
        F_zf=b/(a+b)*m*g;
        F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

        F_zr=a/(a+b)*m*g;
        F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

        F_total=sqrt((Nw*Fx)^2+(F_yr^2));
        F_max=0.7*m*g;

        %F_y derivatives
        dF_yfdphi_yf = F_zf*Dy*By*Cy*cos(Cy*atan(By*phi_yf))/(By^2*phi_yf^2 + 1);
        dF_yrdphi_yr = F_zr*Dy*By*Cy*cos(Cy*atan(By*phi_yr))/(By^2*phi_yr^2 + 1);

        dF_yfdu = dF_yfdphi_yf*dphi_yfda_f*da_fdu;
        dF_yfdv = dF_yfdphi_yf*dphi_yfda_f*da_fdv;
        dF_yfdr = dF_yfdphi_yf*dphi_yfda_f*da_fdr;

        if F_total>F_max
            Fx=F_max/F_total*Fx;

            dF_totaldF_yr = F_yr/F_total;
            dF_yr_dF_yr = F_max/F_total - F_max*F_yr/(dF_totaldF_yr^2);

            F_yr=F_max/F_total*F_yr;

            dF_yrdu = dF_yr_dF_yr*dF_yrdphi_yr*dphi_yrda_r*da_rdu;
            dF_yrdv = dF_yr_dF_yr*dF_yrdphi_yr*dphi_yrda_r*da_rdv;
            dF_yrdr = dF_yr_dF_yr*dF_yrdphi_yr*dphi_yrda_r*da_rdr;
        else
            dF_yrdu = dF_yrdphi_yr*dphi_yrda_r*da_rdu;
            dF_yrdv = dF_yrdphi_yr*dphi_yrda_r*da_rdv;
            dF_yrdr = dF_yrdphi_yr*dphi_yrda_r*da_rdr;
        end

        %X derivatives
        dXdu = cos(psi);
        dXdv = -sin(psi);
        dXdpsi = -u*sin(psi) - v*cos(psi);

        %u derivatives
        dudu = -dF_yfdu*sin(delta)/m;
        dudv = -dF_yfdv*sin(delta)/m + r;
        dudr = -dF_yfdr*sin(delta)/m + v;

        %Y derivatives
        dYdu = sin(psi);
        dYdv = cos(psi);
        dYdpsi = u*cos(psi) - v*sin(psi);

        %v derivatives
        dvdu = dF_yfdu*cos(delta)/m + dF_yrdu/m - r;
        dvdv =dF_yfdv*cos(delta)/m + dF_yrdv/m;
        dvdr = dF_yfdr*cos(delta)/m + dF_yrdr/m - u;

        %psi derivatives
        dpsidr = 1;

        %r derivatives
        drdu = a*dF_yfdu*cos(delta)/Iz - b*dF_yrdu/Iz;
        drdv = a*dF_yfdv*cos(delta)/Iz - b*dF_yrdv/Iz;
        drdr = a*dF_yfdr*cos(delta)/Iz - b*dF_yrdr/Iz;


        %    X          u           Y           v           psi         r
        A = [0          dXdu        0           dXdv        dXdpsi      0;
             0          dudu        0           dudv        0           dudr;
             0          dYdu        0           dYdv        dYdpsi      0;
             0          dvdu        0           dvdv        0           dvdr;
             0          0           0           0           0           dpsidr;
             0          drdu        0           drdv        0           drdr];

    end

    function B = calculate_B(state, input)
        Nw=2;
        f=0.01;
        Iz=2667;
        a=1.35;
        b=1.45;
        By=0.27;
        Cy=1.2;
        Dy=0.7;
        Ey=-1.6;
        Shy=0;
        Svy=0;
        m=1400;
        g=9.806;

        X = state(1); u = state(2);
        Y = state(3); v = state(4);
        psi = state(5); r = state(6);

        if u == 0
           u = 0.001; 
        end
        %generate input functions
        Fx=input(1);
        delta=input(2);

        %slip angle functions in degrees
        a_f=rad2deg(delta-atan2(v+a*r,u));
        a_r=rad2deg(-atan2((v-b*r),u));

        %Alpha derivatives
        da_fddelta = rad2deg(1);


        %Nonlinear Tire Dynamics
        phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
        phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

        %Phi derivatives
        dphi_yfda_f = (1-Ey) + Ey/(By^2*(a_f+Shy)^2 + 1);


        F_zf=b/(a+b)*m*g;
        F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

        F_zr=a/(a+b)*m*g;
        F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

        %F_y derivatives
        dF_yfdphi_yf = F_zf*Dy*By*Cy*cos(Cy*atan(By*phi_yf))/(By^2*phi_yf^2 + 1);

        dF_yfddelta = dF_yfdphi_yf*dphi_yfda_f*da_fddelta;


        % Skidding
        F_total=sqrt((Nw*Fx)^2+(F_yr^2));
        F_max=0.7*m*g;

        if F_total>F_max
            dF_totaldFx  = Fx/F_total;
            dFxdFx = F_max/F_total - F_max*Fx/(dF_totaldFx^2);
        else
            dFxdFx = 1;
        end



        %Linear Aproximation: Approximate a_f and a_r as zero
    %     C_af = F_zf*By*Cy*Dy;
    %     C_ar = F_zr*By*Cy*Dy;
    %     
    %     F_yf = C_af*a_f;
    %     F_yr = C_ar*a_r;

        %F_y derivatives
    %     dF_yfddelta = C_af*da_fddelta;


        %Fx derivatives
        dudFx = Nw*dFxdFx/m;

        %delta derivatives
        duddelta = -F_yf*cos(delta)/m - dF_yfddelta*sin(delta)/m;  
        dvddelta = -F_yf*sin(delta)/m + dF_yfddelta*cos(delta)/m;
        drddelta = -a*F_yf*sin(delta)/Iz + a*dF_yfddelta*cos(delta)/Iz;


        %    Fx     delta
        B = [0      0;
             dudFx  duddelta;
             0      0;
             0      dvddelta;
             0      0;
             0      drddelta];
    end

    function dzdt=bike_odefun(x,u)
        %constants
        Nw=2;
        f=0.01;
        Iz=2667;
        a=1.35;
        b=1.45;
        By=0.27;
        Cy=1.2;
        Dy=0.7;
        Ey=-1.6;
        Shy=0;
        Svy=0;
        m=1400;
        g=9.806;


        %generate input functions
        F_x=u(1);
        delta_f=u(2);

        %slip angle functions in degrees
        a_f=rad2deg(delta_f-atan2(x(4)+a*x(6),x(2)));
        a_r=rad2deg(-atan2((x(4)-b*x(6)),x(2)));

        F_zf=b/(a+b)*m*g;
        F_zr=a/(a+b)*m*g;

        %Linear Aproximation: Approximate a_f and a_r as zero
    %     if linearized
    %         C_af = F_zf*By*Cy*Dy;
    %         C_ar = F_zr*By*Cy*Dy;
    %     
    %         F_yf = C_af*a_f;
    %         F_yr= C_ar*a_r;
    %     else
            %Nonlinear Tire Dynamics
        phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
        phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

        F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;
        F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

        F_total=sqrt((Nw*F_x)^2+(F_yr^2));
        F_max=0.7*m*g;

        if F_total>F_max

            F_x=F_max/F_total*F_x;

            F_yr=F_max/F_total*F_yr;
    %         disp('Oooooooooooooooooooooof')
        end

        %vehicle dynamics
        dzdt= [x(2)*cos(x(5))-x(4)*sin(x(5));...
                  (-f*m*g+Nw*F_x-F_yf*sin(delta_f))/m+x(4)*x(6);...
                  x(2)*sin(x(5))+x(4)*cos(x(5));...
                  (F_yf*cos(delta_f)+F_yr)/m-x(2)*x(6);...
                  x(6);...
                  (F_yf*a*cos(delta_f)-F_yr*b)/Iz];
    end
end