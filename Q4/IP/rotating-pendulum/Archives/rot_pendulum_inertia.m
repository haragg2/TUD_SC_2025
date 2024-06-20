function [dx, y] = rot_pendulum(t, x, u, A, B, C, D, E, F, G, H, varargin)

    g = 9.81;
    theta1 = x(1);
    theta2 = x(3);
    theta1_dot = x(2);
    theta2_dot = x(4);

    % A = b + b1 + b2;
    % B = I1 + m1*c1^2 + m2*l1^2;
    % C = m1*c1 + m2*l1;
    % D = l1*c2*m2;
    % E = b2;
    % F = a;
    % G = m2*c2^2 + I2;
    % H = m2*c2;
    
    M_l = [B, D*cos(theta1 - theta2);
        D*cos(theta1 - theta2), G];
    
    C_l = [A*theta1_dot + D*sin(theta1 - theta2)*theta2_dot^2 - E*theta2_dot; 
           -E*(theta1_dot - theta2_dot) - D*sin(theta1 - theta2)*theta1_dot^2];
    
    G_l = [C*sin(theta1);
            H*sin(theta2)];
    
    B_l = [F*u; 0];
    
    theta_ddot = inv(M_l)*(B_l - C_l - G_l*g);

    % Output equations
    y = [theta1;                         
        theta2];

     % State equations
    dx = [theta1_dot;
          theta_ddot(1,1);
          theta2_dot;
          theta_ddot(2,1)];

end