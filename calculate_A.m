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
    da_fdu = rad2deg((v+a*r)/((v+a*r)^2 + u^2));
    da_fdv = rad2deg(-u/((v+a*r)^2 + u^2));
    da_fdr = rad2deg(-a*u/((v+a*r)^2 + u^2));

    da_rdu = rad2deg((v-b*r)/((v-b*r)^2 + u^2));
    da_rdv = rad2deg(-u/((v-b*r)^2 + u^2));
    da_rdr = rad2deg(b*u/((v-b*r)^2 + u^2));
        
    

%     da_fdu = -(v+a*r)/((v+a*r)^2 + u^2);
%     da_fdv = u/((v+a*r)^2 + u^2);
%     da_fdr = a*u/((v+a*r)^2 + u^2);
%     da_fddelta = 1;
%     
%     da_rdu = (v-b*r)/((v-b*r)^2 + u^2);
%     da_rdv = -u/((v-b*r)^2 + u^2);
%     da_rdr = b*u/((v-b*r)^2 + u^2);
    
%     da_fdu = 0;
%     da_fdv = 0;
%     da_fdr = 0;
%     da_fddelta = 0;
%     
%     da_rdu = 0;
%     da_rdv = 0;
%     da_rdr = 0;
    
    %F_y derivatives
    dF_yfdu = C_af*da_fdu;
    dF_yfdv = C_af*da_fdv;
    dF_yfdr = C_af*da_fdr;
    
    dF_yrdu = C_ar*da_rdu;
    dF_yrdv = C_ar*da_rdv;
    dF_yrdr = C_ar*da_rdr;
    
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