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

    %Nonlinear Tire Dynamics
%     phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
%     phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

    F_zf=b/(a+b)*m*g;
%     F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

    F_zr=a/(a+b)*m*g;
%     F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;
    
    %Linear Aproximation: Approximate a_f and a_r as zero
    C_af = F_zf*By*Cy*Dy;
    C_ar = F_zr*By*Cy*Dy;
    
    F_yf = C_af*a_f;
    F_yr = C_ar*a_r;
    
    %Alpha derivatives
    da_fddelta = rad2deg(1);
    
    %F_y derivatives
    dF_yfddelta = C_af*da_fddelta;
     
    
    %Fx derivatives
    dudFx = Nw/m;
    
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