% Decision variables
z = sdpvar(2, 1);

constraints = [2*z(1) + 5*z(2) + 34 <= 0; z(1) == 7];

[A, b, sense] = convert_constraints(constraints);

% Helper function to convert YALMIP constraints to Gurobi format
function [A, b, sense] = convert_constraints(constraints)
    % This function converts YALMIP constraints to Gurobi format
    % Extract constraints in matrix form
    F = [];
    for i = 1:length(constraints)
        F = [F; getbase(constraints(i))];
    end

    A = F(:, 2:end); % Exclude the first column which is the constant term
    b = -F(:, 1); % The constant term (negated)
    sense = repmat('<', size(A, 1), 1); % All constraints are inequalities

    % Adjust sense for equality constraints
    for i = 1:length(constraints)
        if strcmp(constraints(i).type, '==')
            sense(i) = '=';
        end
    end
end