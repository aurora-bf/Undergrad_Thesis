function [t,y] = rk4(f,a,b,ya,n,len,m)
%len is the length of an initial signal
%m is the level of wavelet decomp
%http://www.math.mcgill.ca/gantumur/math579w10/matlab/rk4.m
%ya is initial condition, n is number of subintervals, a and b are time
%interval
h = (b - a) / n;
halfh = h / 2;
y(1,:) = ya;
t(1) = a;
h6 = h/6;
for i = 1 : n
    t(i+1) = t(i) + h;
    th2 = t(i) + halfh;
    s1 = f(t(i), y(i,:),len,m);
    s2 = f(th2, y(i,:) + halfh * s1,len,m);
    s3 = f(th2, y(i,:) + halfh * s2,len,m);
    s4 = f(t(i+1), y(i,:) + h * s3,len,m);
    y(i+1,:) = y(i,:) + (s1' + s2'+s2' + s3'+s3' + s4') * h6; %i transposed all of those
end
end