function [X, W] = gauss_quadrature_2d(npoints1, npoints2, xia, xib,  etaa, etab)
% Return the gauss quadrauture points and corresponding weights
% https://pomax.github.io/bezierinfo/legendre-gauss.html
%
% npoints1 : number of points in 1-direction
% npoints2 : number of points in 2-direction
% xia, xib = bounds of the element in xi-direction
% etaa, etab = bounds of the element in eta-direction
%

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

%% Get quadrature points in eta dir

if npoints2 == 2
    etahat = [-0.5773502691896257, 0.5773502691896257];
    wetahat = [1.0, 1.0];
elseif npoints2 == 3
    etahat = [0.0000000000000000, -0.7745966692414834, 0.7745966692414834];
    wetahat = [0.8888888888888888,  0.5555555555555556, 0.5555555555555556];
elseif npoints2 == 4
    etahat = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
    wetahat = [0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538];
elseif npoints2 == 5
    etahat = [0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640];
    wetahat = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
else
    error('invalid number of gauss points in eta dir');
end

%% Tranfom points and weights to from etahat to eta domain

wxi = (xib-xia)*wxihat/2.0;
xi  = (xib-xia)*xihat/2.0 + (xib+xia)/2.0;

weta = (etab-etaa)*wetahat/2.0;
eta  = (etab-etaa)*etahat/2.0 + (etab+etaa)/2.0;

X = zeros(npoints1*npoints2,2);
W = zeros(npoints1*npoints2,1);

% Perform tensor product to obtain 2d grid
ctr = 0;
for i = 1 : npoints1
    for j = 1 : npoints2
        ctr = ctr + 1;
        W(ctr)   = wxi(i)*weta(j);
        X(ctr,:) = [xi(i), eta(j)];
    end
end
end
