%
% Implementation of Kalman-like unbiased FIR filter
%
% See page 987 of
% https://www.researchgate.net/publication/233854533_A_Kalman-like_FIR_estimator_ignoring_noise_and_initial_conditions
%
% We're assuming that our A, B, C, and D are time-invariant
%
1; % Don't warn about filename/function-name mismatch

% Set up initial values
%
% UFIR_N is the sliding window size
% A,C from the discrete system state space
function [xhat,F] = ufir_init(A,C,Y,UFIR_N)
    % Let's just try it by initializing all our measurements initially to
    % zero
    %Y = zeros(size(C,1),UFIR_N);

    Csm = zeros(size(C,1)*UFIR_N,size(C,2));
    for i=1:UFIR_N-1
        Csm(size(C,1)*i:size(C,1)*i+size(C,1)-1,:) = C*A^(UFIR_N-1-i);
    end

    Ysm = zeros(size(Y,1)*UFIR_N,1);
    for i=1:UFIR_N-1
        Ysm(size(Y,1)*i:size(Y,1)*i+size(Y,1)-1,1) = Y(:,UFIR_N-i+1);
    end

    K = size(A,1);
    scriptA = A^(K-1);
    P = inv(Csm'*Csm);
    F = scriptA*P*scriptA';
    xhat = scriptA*P*Csm'*Ysm;
end

% Update our state estimates recursively
function [xhat,F] = ufir_update(A,C,y,xhat_last,F_last)
    %upsilon = A; % p = 0, since we're filtering
    F = inv(C'*C + inv(A*F_last*A'));
    %F = pinv(C'*C + pinv(A*F_last*A'));
    %xhat = A*xhat_last + A*inv(upsilon)*F*C'*(y-C*upsilon*xhat_last);
    xhat = A*xhat_last + F*C'*(y-C*A*xhat_last);
end
