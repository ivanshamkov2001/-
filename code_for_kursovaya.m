f_u = @(u,r) exp(r*(1-u.^3))*(1-6*r*u.^3)./(2*sqrt(u));
f = @(u,r) sqrt(u).*exp(r*(1-u.^3));
n = 10;
t = 1:n;
res = zeros(n, 1);
res(1) = 0.01;
r = 0.25;
for i = 2:n
    res(i) = f(res(i-1),r);
end
plot(t,res, 'r*', t,res, 'r')
xlabel('t')
ylabel('u_t')
%legend('$f_u$', 'interpreter', 'latex')
%%
n = 20;
t = 1:n;
res = zeros(n, 1);
res(1) = 1.01;
r = 0.25;
for i = 2:n
    res(i) = f(res(i-1),r);
end
plot(t,res, 'r*', t,res, 'r')
xlabel('t')
ylabel('u_t')
%%
clear all
clc
u = sym('u');
r = sym('r');
f = @(u,r) r*sqrt(u)*(1 - u^3);
df = diff(f)
%%
clear all
clc
u = sym('u');
r = sym('r');
f = @(u,r) r*sqrt(u).*(1-u.^3);
f3 = f(f(f(u,r),r),r);
%plot(linspace(0,2,0.01),f(linspace(0,2,0.01)));
%f3 = f(f(u,r),r);
f3_u = diff(f3, u);
x0 = [2, 3];
%x0 = [0.6, 1.4];
df1u = diff(f3, u); 
df1r = diff(f3, r); 
df2u = diff(f3_u, u); 
df2r = diff(f3_u, r); 
f0 = [double(subs(f3, {u, r}, [x0(1), x0(2)])) - x0(1); double(subs(f3_u, {u, r}, [x0(1), x0(2)])) - 1];
df1u = double(subs(subs(df1u, u, x0(1)), r, x0(2)));
df1r = double(subs(subs(df1r, u, x0(1)), r, x0(2))); 
df2u = double(subs(subs(df2u, u, x0(1)), r, x0(2)));
df2r = double(subs(subs(df2r, u, x0(1)), r, x0(2)));

w = [df1u df1r; df2u df2r];  %матрица Якоби
x = x0 - inv(w) * f0; %решение после первой итерации 
eps = 0.01; 
k = 1;

while (abs(x(1) - x0(1)) > eps) || (abs(x(2) - x0(2)) > eps) 
    x0(1) = x(1);
    x0(2) = x(2);
    f0 = [double(subs(f3, {u, r}, [x0(1), x0(2)])) - x0(1); double(subs(f3_u, {u, r}, [x0(1), x0(2)])) - 1];
    df1u = double(subs(subs(df1u, u, x0(1)), r, x0(2))); 
    df1r = double(subs(subs(df1r, u, x0(1)), r, x0(2))); 
    df2u = double(subs(subs(df2u, u, x0(1)), r, x0(2)));
    df2r = double(subs(subs(df2r, u, x0(1)), r, x0(2)));
    w = [df1u df1r; df2u df2r];
    x = x0 - (inv(w) * f0)';
    k = k + 1; 
    if k > 100
       break;
    end
end
fprintf('r = %f, u = %f\n', x0(2), x0(1));
%% бифуркационная диаграмма
f = @(u,r) sqrt(u).*exp(r*(1-u.^3));
n_r = 400;
n_dots = 100;
res = zeros(n_r * n_dots, 2);
r = 0;
for j = 1:n_r
    r = r + 0.005;
    u0 = 0.1;
    for i = 1:500
        u0 = f(u0, r);
    end
    for i = 1:n_dots
        u0 = f(u0, r);
        res(j * i, 1) = r;
        res(j * i, 2) = u0;
    end
end
hold on
plot(res(:,1), res(:,2), 'r.')
xlabel('r')
ylabel('u')
res = zeros(n_r * n_dots, 2);
r = 0;
for j = 1:n_r
    r = r + 0.005;
    u0 = 3;
    for i = 1:500
        u0 = f(u0, r);
    end
    for i = 1:n_dots
        u0 = f(u0, r);
        res(j * i, 1) = r;
        res(j * i, 2) = u0;
    end
end
hold on
plot(res(:,1), res(:,2), 'r.')
%% Показатель Ляпунова
r = 0;
res = zeros(400, 1);
for i = 1:400
    r = r + 0.005;
    f = @(u) sqrt(u).*exp(r*(1-u.^3));
    f_u = @(u) exp(-r*(u^3 - 1))/(2*u^(1/2)) - 3*r*u^(5/2)*exp(-r*(u^3 - 1));
    u_t = 0.1;
    n = 1000;
    sum = log(abs(f_u(u_t)))/n;
    for j = 1:n
        u_t = f(u_t);
        sum = sum + log(abs(f_u(u_t)))/n;
    end
    res(i) = sum;
end
plot(linspace(0.005, 2, 400), res, 'r', [0, 2], [0 0], 'k')
xlabel('r')
ylabel('h')
%% устойчивость (1,1) при r < 1/48
f = @(u,v,r) sqrt(u).*exp(r*(1-v.^3));
n = 50;
t = 1:n;
res = zeros(n, 1);
res(1) = 1.1;
res(2) = res(1);
r = 1/96;
for i = 3:n
    res(i) = f(res(i-1),res(i-2),r);
end
plot(t,res, 'r*', t,res, 'r')
xlabel('t')
ylabel('u_t')
%% устойчивость (1,1)
f = @(u,v,r) sqrt(u).*exp(r*(1-v.^3));
n = 50;
t = 1:n;
res = zeros(n, 1);
res(1) = 1.1;
res(2) = res(1);
r = 0.25;
for i = 3:n
    res(i) = f(res(i-1),res(i-2),r);
end
plot(t,res, 'r*', t,res, 'r')
xlabel('t')
ylabel('u_t')
%% неустойчивость (1,1)
f = @(u,v,r) sqrt(u).*exp(r*(1-v.^3));
n = 50;
t = 1:n;
res = zeros(n, 1);
res(1) = 1.1;
res(2) = res(1);
r = 0.5;
for i = 3:n
    res(i) = f(res(i-1),res(i-2),r);
end
plot(t,res, 'r*', t,res, 'r')
xlabel('t')
ylabel('u_t')
%% Неймарк - Сакер
r = 1/3;
f = @(u,v,r) sqrt(u).*exp(r*(1-v.^3));
u0 = 1.1;
v0 = u0;
n = 5000;
res = zeros(n,2);
for i = 1:n
    t = f(u0, v0, r);
    v0 = u0;
    u0 = t;
    res(i,1) = u0;
    res(i,2) = v0;
end
plot(res(:,1), res(:,2), 'r.')
xlabel('u')
ylabel('v')
%%
r = 1/3 - 0.001;
f = @(u,v,r) sqrt(u).*exp(r*(1-v.^3));
u0 = 1.1;
v0 = u0;
n = 5000;
res = zeros(n,2);
for i = 1:n
    t = f(u0, v0, r);
    v0 = u0;
    u0 = t;
    res(i,1) = u0;
    res(i,2) = v0;
end
plot(res(:,1), res(:,2), 'r.')
xlabel('u')
ylabel('v')
%%
r = 1/3 + 0.001;
f = @(u,v,r) sqrt(u).*exp(r*(1-v.^3));
u0 = 1.01;
v0 = u0;
n = 5000;
res = zeros(n,2);
for i = 1:n
    t = f(u0, v0, r);
    v0 = u0;
    u0 = t;
    res(i,1) = u0;
    res(i,2) = v0;
end
plot(res(:,1), res(:,2), 'r.')
xlabel('u')
ylabel('v')
%%
r = 1/3 + 0.001;
f = @(u,v,r) sqrt(u).*exp(r*(1-v.^3));
u0 = 1.1;
v0 = u0;
n = 5000;
res = zeros(n,2);
for i = 1:n
    t = f(u0, v0, r);
    v0 = u0;
    u0 = t;
    res(i,1) = u0;
    res(i,2) = v0;
end
plot(res(:,1), res(:,2), 'r.')
xlabel('u')
ylabel('v')