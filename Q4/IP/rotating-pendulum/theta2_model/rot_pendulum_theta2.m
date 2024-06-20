function [dx, y] = rot_pendulum(t, x, u, m2, c2, b2, I2, varargin)

    g = 9.81;
    theta2 = x(1);
    theta2_dot = x(2);

    theta2_ddot = (-m2*c2/(I2+m2*c2^2))*g*sin(theta2) - (b2/(I2+m2*c2^2))*theta2_dot;

    % Output equations
    y = [theta2];

     % State equations
    dx = [theta2_dot;
          theta2_ddot];

end