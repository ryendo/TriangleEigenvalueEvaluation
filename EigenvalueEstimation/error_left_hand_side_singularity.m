function val = error_left_hand_side_singularity(sigma,kappa)
    
    kappa = kappa(:)';

    abs_error = (I_intval(I_sup(sigma))-I_intval(I_inf(sigma)))*intval(sup(abs(partialF_kappa_sigma_singularity(kappa, sigma))))./intval(inf(abs(partialF_kappa_kappa_singularity(kappa, sigma))));
    val = kappa + I_hull(-abs_error,abs_error);
    
    val = val(:)';
end


function val = partialF_kappa_kappa_singularity(kappa, sigma)
    % partialF_kappa_s_new calculates the partial derivative of F with respect to \kappa (vectorized).
    %
    % Input:
    % kappa - parameter \kappa (vector)
    % s - parameter s (vector or scalar)
    %
    % Output:
    % val - value of \partial F / \partial \kappa (vector)

    % Ensure kappa and s are column vectors
    kappa = kappa(:)'; 
    sigma = sigma(:)'; % s transposed to support broadcasting

    gi = isa(kappa, 'gradient');
    c1 = GI_intval(1, gi);
    c2 = GI_intval(2, gi);
    c3 = GI_intval(3, gi);

    % Calculate u_+ and u_-
    u_plus = (c2 - sigma^c3).^(c2 / c3) .* kappa;
    u_minus = sigma.^c2 .* kappa;

    % Modified Airy function evaluations
    A_plus = I_minus_Ai(u_plus);
    A_plus_p = I_d_minus_Ai(u_plus);
    A_minus = I_minus_Ai(u_minus);
    A_minus_p = I_d_minus_Ai(u_minus);

    % Compute the terms of the partial derivative
    term1 = A_plus_p .* A_minus_p;
    term2 = (kappa .*(c2 - sigma.^c3).^(c1 / c3)).* sigma .*  A_plus .* A_minus;

    % Combine terms
    val = c2 .*  (term1 + term2);
end

function val = partialF_kappa_sigma_singularity(kappa, sigma)
    % partialF_kappa_s_new calculates the partial derivative of F with respect to \kappa (vectorized).
    %
    % Input:
    % kappa - parameter \kappa (vector)
    % s - parameter s (vector or scalar)
    %
    % Output:
    % val - value of \partial F / \partial \kappa (vector)

    % Ensure kappa and s are column vectors
    kappa = kappa(:)'; 
    sigma = sigma(:)'; % s transposed to support broadcasting

    gi = isa(kappa, 'gradient');
    c1 = GI_intval(1, gi);
    c2 = GI_intval(2, gi);
    c3 = GI_intval(3, gi);
    c4 = GI_intval(4, gi);

    % Calculate u_+ and u_-
    u_plus = (c2 - sigma^c3).^(c2 / c3) .* kappa;
    u_minus = sigma.^c2 .* kappa;

    % Modified Airy function evaluations
    A_plus = I_minus_Ai(u_plus);
    A_plus_p = I_d_minus_Ai(u_plus);
    A_minus = I_minus_Ai(u_minus);
    A_minus_p = I_d_minus_Ai(u_minus);

    % % Compute the terms of the partial derivative
    % term1 = -sigma^2 .*(c2 - sigma^3).^(-c2 / c3) .* A_plus .* A_minus_p;
    % term2 = (c2 - sigma^3).^(c1 / c3) .* A_minus .* A_plus_p;
    % term3 = c4 .* kappa .* A_plus_p .* A_minus_p;
    % term4 = c4 .* ((c1 - sigma.^c2).^(c1 / c3)) .* kappa.^c2 .* A_plus .* A_minus;
    % 
    % % Combine terms
    % val = (c1 / c3) .* (term1 + term2 + term3 + term4);

    val = -sigma.^c2 .* (c2 - sigma.^c3).^(- c2 / c3) .* A_plus .* A_minus_p ...
    + (c2 - sigma.^c3).^(c1 / c3) .* (-2 * sigma.^c2 .* (c2 - sigma.^c3).^(-c1/ c3) .* kappa .* A_plus_p .* A_minus_p ...
    + 2 * sigma.^c3 .* kappa.^2 .* A_plus .* A_minus)...
    + A_minus .* A_plus_p + c2.* kappa.* sigma.^2 .* A_plus_p .* A_minus_p...
    + -2 * sigma.^c3 .* (c2 - sigma.^c3).^(c1/ c3) .* kappa.^2 .* A_plus .* A_minus;
end
