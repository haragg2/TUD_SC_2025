function [th_ekf,P_ekf] = measupdate_ekf(th_k,P_k,y_k,mic_locations,R)
    % INPUT
    % th_k              mean of predictive distribution for k|k-1
    % P_k               covariance of predictive distribution for k|k-1
    % y_k               kth measurment
    % mic_locations     microphone locations
    % R                 measurement noise covariance matrix
    % OUTPUT
    % th_ekf            mean of filtering distribution for k|k
    % p_ekf             covariance matrix of filtering distribution for k|k

    c               = 343; %In meters per second

end