clear all
close all
clc


tic
nt = 60; %Tempo total em (s)
dt = 0.00001; %intervalo de tempo para convergencia (s)


L = 0.03; %comprimento total (m)
nx = 100; %n�mero de n�s em x
dx = L/nx; % comprimento da c�lula em x

k = 27; %condutividade termica (W/mK)
ro = 7805; %densidade (kg/m3)
c = 600; %calor especifico (J/kgK)
alfa = k/(ro*c); %difusidade t�rmica (m�/s)
To = 30; %Temperatura inicial (C�)
fluxo = 1e6; %fluxo de calor (W/m�)
T = ones(nx,1)*To; %Vetor de temperatura inicial (C�)

m=0; %Varia��o auxiliar na programa��o
tolerancia = 0.01; %Crit�rio de toler�ncia para convergencia da malha

%% L�gica cria��o arquivo .txt
str1 = 'T';
str2 = num2str(nt);
str3 = 's.txt';
filename = strcat(str1,str2,str3);
arq = fopen(filename, 'wt');
%%

for t=dt:dt:nt
    
    T_old = T; %Criando vetor temperatura Tp^n
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% C�lula 1: x=dx/2
    m = m+1;
    Tw = (fluxo*dx/k + T_old(m));
    %Calculo do Tp^(n+1)
    T(m) = T_old(m) + alfa*dt*((Tw-2*T_old(m)+T_old(m+1))/(dx^2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% C�lulas Centrais: dx/2 < x > L-dx/2
    for i=2:nx-1
        m = m+1;
        %Calculo Tp^(n+1)
        T(m) = T_old(m) + alfa*dt*((T_old(m-1)-2*T_old(m)+T_old(m+1))/(dx^2));
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% C�lula 3 = x=L-dx/2
    m = m+1;
    Te = T_old(m);
    %Calculando Tp^(n+1)
    Tm = T_old(m) + alfa*dt*((T_old(m-1)-2*T_old(m)+Te)/(dx^2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    m=0;
    
end


x = dx/2:dx:L-dx/2

%imprimir os dados em .txt
%for i=1:nx
%%end

toc

plot(x,T)
title('Gr�fico de temperatura em 60s')
xlabel('Comprimento (m)')
ylabel('Temperatura (C�)')
grid


     
    
