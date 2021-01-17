function [x,d,r] = molecule_NCGM(x,k,sigma,epsilon,kz)
%molecule_NCGM performs the Nonlinear Conjugate Gradient Method to find the
%    minimum of the energy of an N particle system.
%
% Arguments: x (3Nx1 real vector) gives the positions of the atoms.
%            k (real constant) gives the harmonic potential constant.
%            sigma (real constant) gives the sigma coulomb constant.
%            epsilon (real constant) gives the epsilon coulomb constant.
%            kz (real constant) gives the extra term for z harmonic.
%
% Returns: x (3Nx1 real vector) gives the positions of the atoms at minimum
%          d (3Nx1 real vector) gives the d vector minimum
%          r (3Nx1 real vector) gives the r vector at minimum
tic
e1 = 1e-6;
e2 = 1e-12;
r = -gradE(x,k,sigma,epsilon,kz);
d = r;
while norm(r) > e1
    alpha = 1/norm(d);
    while norm(alpha*d) > e2
        df2d = f2d(x,k,sigma,epsilon,kz,d);
        if df2d > 0
            alpha = -gradE(x,k,sigma,epsilon,kz)'*d/df2d;
            x = x + alpha*d;
        else
            e1 = e1*2; % larger tolerance to prevent further df2d<0 or NaN
            e2 = e2*2;
            x = x - e1*r;
            alpha = 0; % disp('beware: df2d < 0 or NaN')
        end
    end
    r1 = -gradE(x,k,sigma,epsilon,kz);
    beta = r1'*r1/(r'*r);
    d = r1 + beta*d;
    r = r1;
end
toc
end

