function [commonFactor, normalizedK, normalizedM] = normalizeMatrices(K, M)
    % Normalize two matrices K and M by dividing both with a common factor
    % The common factor minimizes the absolute value of powers of symbolic variables t and s
    % The equivalence of the generalized eigenvalue problem Kx = lambda Mx is preserved

    syms t s real

    % Flatten the elements of both matrices into vectors
    elementsK = reshape(K, [], 1);
    elementsM = reshape(M, [], 1);

    % Combine all elements from both matrices into a single vector
    allElements = [elementsK; elementsM];

    % Compute the degree of each element with respect to t and s
    exponentsT = abs(feval(symengine, 'degree', allElements, t)); % Degrees with respect to t
    exponentsS = abs(feval(symengine, 'degree', allElements, s)); % Degrees with respect to s

    % Sum the absolute values of the degrees of t and s
    totalExponents = exponentsT + exponentsS;

    % Find the element with the smallest total degree
    [~, idx] = min(totalExponents);
    commonFactor = allElements(idx); % Select the common factor

    % Normalize the matrices by dividing by the common factor
    normalizedK = simplify(K / commonFactor);
    normalizedM = simplify(M / commonFactor);
end
