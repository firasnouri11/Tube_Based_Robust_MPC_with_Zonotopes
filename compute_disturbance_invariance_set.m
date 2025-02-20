function [Error_Zonotope, Xc_robust, Uc_robust] = compute_disturbance_invariance_set(A,B,K,W,s,Xc,Uc)
    % calculate the bounding set Epsilon as an outer approximation of the
    % MRPI set satisfying minkowski sum of A_cl*Epsilon and disturbance W as a subset of Epsilon 
    A_cl = A - B .* K;
    Epsilon = W;
    for j = 1:s
        Epsilon = Epsilon + A_cl^j * W;
    end
    Error_Zonotope = Epsilon;
    Xc_robust = minkDiff(Xc,Epsilon);
    Uc_robust = minkDiff(Uc,K * Epsilon);
end