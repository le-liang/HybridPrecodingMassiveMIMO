% expintn exponential intergral of the order n
% r = expintn(x, n) gives integral from 1 to inf of (exp(-xt)/t^n)dt, for x > 0
% Reference: 
%      I. S. Gradshteyn and I. M. Ryzhik, Table of Integrals, Series, and Products
%      Eq. (8.352.5)
% By Le Liang, UVic, Dec. 05, 2013

function r = expintn_noMaple(x, n)

if n >= 2
    sum = 0;
    for m = 0 : (n-2)
        sum = sum + (-1)^m*factorial(m)/x^(m+1);
    end

    minuend = x^(n-1) * (-1)^(n+1)/factorial(n-1) * expint(x);
    subtractor = x^(n-1) * (-1)^(n+1)/factorial(n-1) * exp(-x) * sum;
    r = minuend - subtractor;
else
    r = expint(x);
end