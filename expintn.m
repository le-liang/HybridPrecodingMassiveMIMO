% expintn(x, n) is an an exponential integral of order n with input x,
% computed using maple

% By Le Liang, UVic

function r = expintn(x, n)

cmd_str = sprintf('evalf[16](Ei(%d, %g));', n, x);
r = double(maple(cmd_str));
