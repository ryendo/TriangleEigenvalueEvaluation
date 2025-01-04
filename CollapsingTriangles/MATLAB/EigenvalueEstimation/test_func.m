function val = test_func(s)
        % partialF_kappa_s calculates the partial derivative of F with respect to s.
    %
    % Input:
    % kappa - parameter \kappa (vector or scalar)
    % s - parameter s (vector or scalar, must be broadcastable with kappa)
    %
    % Output:
    % val - value of \partial F / \partial s (vector or scalar)

    % Ensure kappa and s are row vectors for element-wise operations
    kappa = [1,2];
    kappa = kappa(:)';
    s = s(:)';

    % Constants
    c1 = I_intval(1);
    c2 = I_intval(2);
    c3 = I_intval(3);
    c4 = I_intval(4);

    % Calculate common terms
    factor1 = (c1 + s).^(c2 / c3);
    factor2 = (c1 - s).^(c2 / c3);
    factor3 = (c1 - s.^c2).^(c1 / c3);

    % Modified Airy function evaluations
    A1 = I_minus_Ai(factor1 .* kappa);
    A1p = I_d_minus_Ai(factor1 .* kappa);
    A2 = I_minus_Ai(factor2 .* kappa);
    A2p = I_d_minus_Ai(factor2 .* kappa);

    % Compute the terms of the partial derivative
    term1 = (c1 + s).^(-c2 / c3) .* A1 .* A2p;
    term2 = (c1 - s).^(-c2 / c3) .* A2 .* A1p;
    term3 = c4 .* kappa .* A1p .* A2p;
    term4 = c4 .* factor3 .* kappa.^c2 .* A1 .* A2;

    % Combine terms
    val = (c1 / c3) .* (term1 + term2 + term3 + term4);

    % Ensure val matches the shape of the inputs
    val = reshape(val, size(kappa));
end

