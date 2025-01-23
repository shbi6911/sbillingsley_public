


function [Kp, Kd]= timeConst2gains(tau1, tau2, I)
        lambda1=-1/tau1;
        lambda2=-1/tau2;
    % Creating linear system Ax=b, where x= [Kp; Kd]
        A= [1 lambda1; ...
            1 lambda2] ;
        b= -I*[ lambda1^2; lambda2^2 ]; % solution vector
    % Inverting and LH multiplying A to b
        x= (A^-1)*b;
    % Proportional gain, Kp
        Kp= x(1);
    % Derivative gain, Kd
        Kd= x(2);
end