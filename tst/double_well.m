clear all;
clc;
%% initialize parameters

par.a   = 1.0;
par.b   = 1.0;
par.K   = 0.0;
par.D0  = 0.1;
par.tau = 1.0;

par.init.x0 = 0.9*sqrt(2*par.b/par.a);
par.init.s  = 0.1;

par.dim.t   = 1;
par.dim.tau = 1;
par.dim.all = par.dim.t+par.dim.tau;

%% initialize Xt grid

X{1}.min = -3;
X{1}.max =  3;
X{1}.pts = 2^9;

X{1}.len = X{1}.max-X{1}.min;
X{1}.del = X{1}.len/(X{1}.pts-1);

X{1}.S = linspace(X{1}.min,X{1}.max,X{1}.pts)';
X{1}.G = X{1}.S;

%% initialize Xtau grid

Y{1}.min = -3;
Y{1}.max =  3;
Y{1}.pts = 2^4;

Y{1}.len = Y{1}.max-Y{1}.min;
Y{1}.del = Y{1}.len/(Y{1}.pts-1);

Y{1}.S = linspace(Y{1}.min,Y{1}.max,Y{1}.pts)';
Y{1}.G = Y{1}.S;

%% initialize compound X grid

Z{1} = X{1};
Z{2} = Y{1};

[Z{1}.G,Z{2}.G] = meshgrid(X{1}.G,Y{1}.G);

%% initialize time

T.min =  0.0;
T.max = 10.0;
T.pts = 20;

T.len = T.max-T.min;
T.del = T.len/(T.pts-1);

T.S = linspace(T.min,T.max,T.pts)';

T.hist.pts = round(par.tau/T.del);
par.actual_tau = T.hist.pts*T.del;

T.hist.idx = 1;

%% initialize rho 1

P{1}.G = init_gauss(X{1}.G,par.init);

P{1}.hist = cell(T.hist.pts,1);

for itau = 1:T.hist.pts
    P{1}.hist{itau} = P{1}.G;
end

%% initialize rho 2

P{2}.G = zeros(size(Z{1}.G));

%{
P{1}.hist = cell(T.hist.pts,1);

for itau = 1:T.hist.pts
    P{1}.hist{itau} = P{1}.G;
end
%}

%% 




%% visual check potential

V = potential(X{1}.G,par);

plot(X{1}.G,P{1}.G);
hold on;
plot(X{1}.G,V);
hold off;

%% calculate a rhs for rho1 at t=t_i

x = X{1}.G;
p = par;

rhs1 = -force1_dx(x,p) .* P{1}.G
%rhs2 = -force1(x,p) .* dP???



%%

function V = potential(X,par)
    a = par.a;
    b = par.b;
    V = a/4*X.^4-b/2*X.^2;
end

function F1 = force1(X,par)
    a = par.a;
    b = par.b;
    K = par.K;
    F1 = -a*X.^3 + (b+K)*X;
end

function F2 = force2(X,par)
    K = par.K;
    F2 = -K*X;
end

function F1dx = force1_dx(X,par)
    a = par.a;
    b = par.b;
    K = par.K;
    F1dx = -3*a*X.^2 + (b+K);
end


function P0 = init_gauss(X,init_par)
    x0 = init_par.x0;
    s = init_par.s;

    P0 = 1/(s*sqrt(2*pi))*exp(-1/(2*s^2)*(X-x0).^2);
end
