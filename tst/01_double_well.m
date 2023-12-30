clear all;
clc;
%%

X{1}.min = -5;
X{1}.max =  5;
X{1}.pts = 2^5;

X{1}.len = X{1}.max-X{1}.min;
X{1}.del = X{1}.len/(X{1}.pts-1);

X{1}.S = linspace(X{1}.min,X{1}.max,X{1}.pts)';
X{1}.G = X{1}.S;

%%

Y{1}.min = -5;
Y{1}.max =  5;
Y{1}.pts = 2^5;

Y{1}.len = Y{1}.max-Y{1}.min;
Y{1}.del = Y{1}.len/(Y{1}.pts-1);

Y{1}.S = linspace(Y{1}.min,Y{1}.max,Y{1}.pts)';
Y{1}.G = Y{1}.S;

%%

Z{1} = X{1};
Z{2} = Y{1};

Z{1}.G = 
