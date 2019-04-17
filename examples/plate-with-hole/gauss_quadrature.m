function [X, W] = gauss_quadrature(npoints)
% Return the gauss quadrauture points and corresponding weights
% https://pomax.github.io/bezierinfo/legendre-gauss.html
if npoints == 2
    x = [-0.5773502691896257, 0.5773502691896257];
    w = [1.0, 1.0];
elseif npoints == 3
    x = [0.0000000000000000, -0.7745966692414834, 0.7745966692414834];
    w = [0.8888888888888888,  0.5555555555555556, 0.5555555555555556];
elseif npoints == 4
    x = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
    w = [0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538];
elseif npoints == 5
    x = [0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640];
    w = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
else
    error('invalid number of gauss points');
end

X = zeros(npoints,2);
W = zeros(npoints,1);

% Perform tensor product to obtain 2d grid
ctr = 0;
for i = 1 : npoints
    for j = 1 : npoints
        ctr = ctr + 1;
        W(ctr)   = w(i)*w(j);
        X(ctr,:) = [x(i), x(j)];        
    end
end
end
