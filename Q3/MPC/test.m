m_w = 1.5;          % kg; mass of wheel
m_p = 6.16;         % kg; mass of pendulum
l = 0.46;           % m; distance between COG of pendulum and center of wheel
R = 0.054;          % m; radius of wheel
g = 9.81;           % m/s^2; gravity
J = 0.1;            % kgm^2; wheel centroidal inertia?
 
x = 1.452;
a = J + ((m_w + m_p) * R^2);
c =  m_p * l^2;
d = (m_w + m_p) * g * R * sin(beta);

compute_expression(x, 0, 10*pi/180, 0, a, c, l, m_p, R, g, d)

function result = compute_expression(x, x_dot, beta, tau, a, c, l, m, r, g, d)
    numerator_1 = (c * l * m * r * cos(x) * x_dot^2 + g * l^2 * m^2 * r * sin(x) * sin(beta + x) ...
        - g * l^2 * m^2 * r * cos(x) * cos(beta + x) - l * m * r * tau * sin(beta + x));
    denominator_1 = (a * c - l^2 * m^2 * r^2 * cos(beta + x)^2);
    
    numerator_2 = -2 * l^2 * m^2 * r^2 * sin(beta + x) * cos(beta + x) * (-c * d + c * l * m * r * sin(x) * x_dot^2 ...
        + tau * (c + l * m * r * cos(beta + x)) - g * l^2 * m^2 * r * sin(x) * cos(beta + x));
    denominator_2 = (a * c - l^2 * m^2 * r^2 * cos(beta + x)^2)^2;
    
    result = (numerator_1 / denominator_1) - (numerator_2 / denominator_2);
end
