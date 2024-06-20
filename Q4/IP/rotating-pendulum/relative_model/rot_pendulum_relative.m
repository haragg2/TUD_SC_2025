function [dx, y] = rot_pendulum(t, x, u, A, C, D, F, a, varargin)
    g = 9.81;
    
    % After SYS ID of Theta 2
    I2 = 0.000110802;
    m2 = 0.05914;
    c2 = 0.0563237;
    b2 = 5.65455e-05;

    theta1 = x(1);
    theta2 = x(3);
    theta1_dot = x(2);
    theta2_dot = x(4);

    %A = I1 + m1*c1^2
    %C = l1
    %D = m1*c1 
    B = I2 + m2*c2^2;
    E = m2*c2;
    % F = b + b1;
    
    M_l = [A + m2*C^2 + B + 2*m2*C*c2*cos(theta2), B + m2*C*c2*cos(theta2);
           B + m2*C*c2*cos(theta2), B];
    
    C_l = [F - 2*m2*C*c2*theta2_dot*sin(theta2), -m2*C*c2*theta2_dot*sin(theta2); 
           m2*C*c2*theta1_dot*sin(theta2), b2];
    
    G_l = [(D + m2*C)*sin(theta1) + E*sin(theta1+theta2);
           E*sin(theta1+theta2)];

    B_l = [a*u; 0];
    
    theta_ddot = inv(M_l)*(B_l - C_l*[theta1_dot; theta2_dot] - G_l*g);

    % Output equations
    y = [theta1;                         
        theta2];

     % State equations
    dx = [theta1_dot;
          theta_ddot(1,1);
          theta2_dot;
          theta_ddot(2,1)];

end