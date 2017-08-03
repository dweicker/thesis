n = 9;
x = [-1.0, -0.8997579954114601573123, -0.6771862795107377534459, -0.3631174638261781587108, 0.0, 0.3631174638261781587108, 0.6771862795107377534459, 0.8997579954114601573123, 1.0];
w = [0.02777777777777777777778, 0.165495361560805525046, 0.2745387125001617352807, 0.346428510973046345115, 0.3715192743764172335601, 0.346428510973046345115, 0.2745387125001617352807, 0.165495361560805525046, 0.02777777777777777777778];
L = [126.026890 -81.542936 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ;
-242.909599 438.000000 -242.909599 30.168178 -8.706415 3.657143 -1.900603 1.117610 -0.676313 0.250000 0.000000; 
0.000000 -81.542936 126.026890 -51.999472 10.047938 -3.701488 1.814357 -1.035872 0.617616 -0.227033 0.000000 ;
0.000000 6.104822 -31.345930 44.327981 -22.776416 5.073441 -2.075939 1.090319 -0.624437 0.226159 0.000000 ;
0.000000 -1.396218 4.800087 -18.049923 27.645129 -15.707953 3.792065 -1.645147 0.866752 -0.304793 0.000000 ;
0.000000 0.546875 -1.648849 3.749081 -14.647107 24.000000 -14.647107 3.749081 -1.648849 0.546875 0.000000 ;
0.000000 -0.304793 0.866752 -1.645147 3.792065 -15.707953 27.645129 -18.049923 4.800087 -1.396218 0.000000 ;
0.000000 0.226159 -0.624437 1.090319 -2.075939 5.073441 -22.776416 44.327981 -31.345930 6.104822 0.000000 ;
0.000000 -0.227033 0.617616 -1.035872 1.814357 -3.701488 10.047938 -51.999472 126.026890 -81.542936 0.000000 ;
0.000000 0.250000 -0.676313 1.117610 -1.900603 3.657143 -8.706415 30.168178 -242.909599 438.000000 -242.909599 ;
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 -81.542936 126.026890 ];
M = [0.1655    0.0556    0.1655    0.2745    0.3464    0.3715    0.3464    0.2745    0.1655    0.0556 0.1655];

e = ones(n,1);
xx = kron(x',e');
yy = kron(e,x);
rr = zeros(n+2,n+2);
for i=2:n+1
    for j=2:n+1
        rr(i,j) = w(i-1)*w(j-1)*pi*pi*cos(pi*x(i-1)/2)*cos(pi*x(j-1)/2)/(2*M(i)*M(j));
    end
end
rr(2:n+1,1) = 2*rr(2:n+1,2)-rr(2:n+1,3);
rr(1,2:n+1) = 2*rr(2,2:n+1)-rr(3,2:n+1);
rr(2:n+1,n+2)= 2*rr(2:n+1,n+1)-rr(2:n+1,n);
rr(n+2,2:n+1) = 2*rr(n+1,2:n+1)-rr(n,2:n+1);

[V_inv,lambda] = eig(L);
lambda = diag(lambda);
D = zeros(n+2);
for i=1:n+2
    for j = 1:n+2
        D(i,j) = 1/(lambda(i)+lambda(j));
    end
end
V = inv(V_inv);
VRV = V*rr*V_inv;
W = D.*VRV;
z = V_inv*W*V;


