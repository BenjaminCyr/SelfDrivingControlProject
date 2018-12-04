function [A, B] = calculate_jacobian(state, input)
    A = calculate_A(state, input);
    B = calculate_B(state, input);
end