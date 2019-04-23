function [xi, wxi] = gauss_quadrature_1d(npoints1, xia, xib)
% Return the gauss quadrauture points and corresponding weights
% https://pomax.github.io/bezierinfo/legendre-gauss.html
%
% npoints1 : number of points in 1-direction
% xia, xib = bounds of the element in xi-direction

%% Get quadrature points in xi dir
if npoints1 == 2
    xihat = [-0.5773502691896257, 0.5773502691896257];
    wxihat = [1.0, 1.0];
elseif npoints1 == 3
    xihat = [0.0000000000000000, -0.7745966692414834, 0.7745966692414834];
    wxihat = [0.8888888888888888,  0.5555555555555556, 0.5555555555555556];
elseif npoints1 == 4
    xihat = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
    wxihat = [0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538];
elseif npoints1 == 5
    xihat = [0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640];
    wxihat = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
else
    error('invalid number of gauss points in xi dir');
end

%% Tranfom points and weights to from etahat to eta domain
wxi = (xib-xia)*wxihat/2.0;
xi  = (xib-xia)*xihat/2.0 + (xib+xia)/2.0;

end
