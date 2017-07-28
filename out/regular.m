% This solve a small linear system with row ordering

n = 9;
h = 2/(n-1);
z = linspace(-1,1,n)';

e = ones(n,1);
ee = ones(n*n,1);
x = kron(e,z);
y = kron(z,e);

b = h*h*pi*pi*cos(pi*x/2).*cos(pi*y/2)/2;
A = spdiags([-1*ee -1*ee 4*ee -1*ee -1*ee],[-n -1 0 1 n],n*n,n*n);
A(1:n,1:n*n) = speye(n,n*n); b(1:n) = zeros(n,1);
left = (1:n-2)*n + 1;
A(left,1:n*n) = zeros(n-2,n*n); b(left) = zeros(n-2,1);
A(left,left) = speye(n-2);
right = (2:n-1)*n;
A(right,1:n*n) = zeros(n-2,n*n); b(right) = zeros(n-2,1);
A(right,right) = speye(n-2);
A(n*n-n+1:n*n,1:n*n) = zeros(n,n*n); b(n*n-n+1:n*n) = zeros(n,1);
A(n*n-n+1:n*n,n*n-n+1:n*n) = speye(n);


u = A\b;