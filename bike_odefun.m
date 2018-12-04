function dzdt=bike_odefun(x,u, linearized)
    if nargin < 3
       linearized = true; 
    end
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
    if linearized
        C_af = F_zf*By*Cy*Dy;
        C_ar = F_zr*By*Cy*Dy;
    
        F_yf = C_af*a_f;
        F_yr= C_ar*a_r;
    else
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
        end
    end

    %vehicle dynamics
    dzdt= [x(2)*cos(x(5))-x(4)*sin(x(5));...
              (-f*m*g+Nw*F_x-F_yf*sin(delta_f))/m+x(4)*x(6);...
              x(2)*sin(x(5))+x(4)*cos(x(5));...
              (F_yf*cos(delta_f)+F_yr)/m-x(2)*x(6);...
              x(6);...
              (F_yf*a*cos(delta_f)-F_yr*b)/Iz];
end