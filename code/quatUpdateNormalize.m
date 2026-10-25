function q = quatUpdateNormalize(x, q_hat)
    dq = [1 x(1:3, :)'];
    dq = dq / norm(dq);
    q = quatmultiply(q_hat', dq);
    q  = q / norm(q);
end