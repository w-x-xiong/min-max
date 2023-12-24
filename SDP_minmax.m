function [x_vec] = SDP_minmax(Tx, Rx, Rg, rho)

[~, M] = size(Tx);
[~, L] = size(Rx);


cvx_begin

variables dmt_vec(M) gmt_vec(M) gls_vec(L) x_vec(2) z_vec(2) gamma_vec(2) t lambda

minimize t

subject to

z_vec == x_vec + gamma_vec;

[t*eye(2), gamma_vec; gamma_vec', t] == semidefinite(3);

t >= 0;

for m = 1:M
    for l = 1:L
        (Rg(m,l) - rho)^2 + gmt_vec(m) + 2*(Rg(m,l) - rho)*dmt_vec(m) <= gls_vec(l);
        (Rg(m,l) + rho)^2 + gmt_vec(m) + 2*(Rg(m,l) + rho)*dmt_vec(m) >= gls_vec(l);
    end
end

for m = 1:M
    gmt_vec(m) == [Tx(:,m);-1]'*[eye(2),z_vec;z_vec',lambda]*[Tx(:,m);-1];
    [1, dmt_vec(m);dmt_vec(m), gmt_vec(m)] == semidefinite(2);
end

for l = 1:L
    gls_vec(l) == [Rx(:,l);-1]'*[eye(2),z_vec;z_vec',lambda]*[Rx(:,l);-1];
end

[eye(2), z_vec; z_vec', lambda] == semidefinite(3);

cvx_end

end

