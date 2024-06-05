function [dx, y] = rot_pendulum(t, x, u, m2, c2, b2, I2, varargin)
    % I1 = 0.074;
    % I2 = 0.00012;
    % m1 = 0.125;
    % m2 = 0.05;
    % c1 = 0.04; % could be negative
    % c2 = 0.06;
    % l1 = 0.1;
    g = 9.81;
    % b1 = 4.8;
    % b2 = 0.0002;
    
    theta2 = x(1);
    theta2_dot = x(2);

    % A = b + b1 + b2;
    % B = I1 + m1*c1^2 + m2*l1^2;
    % C = m1*c1 + m2*l1;
    % D = l1*c2*m2;
    % E = b2;
    % F = a;
    % G = m2*c2^2 + I2;
    % H = m2*c2;
    
    theta2_ddot = (-m2*c2/(I2+m2*c2^2))*g*sin(theta2) - (b2/(I2+m2*c2^2))*theta2_dot;

    % A = -+m2*c2/(I2+m2*c2^2)
    % B = -b2/(I2+m2*c2^2)

    % Output equations.
    y = [theta2];

     % State equations.
    dx = [theta2_dot;
          theta2_ddot];

end