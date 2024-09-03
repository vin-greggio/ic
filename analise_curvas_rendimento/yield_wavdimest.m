
logFerRate = csvread('parte1.csv');

%Esses códigos foram baseados nos códigos fornecidos por Rodney Fonseca


n = 57;
nt = 64;
u = linspace(0,5.886104031450156,nt);  % nt valores igualmente espaçados de 0 a 5.88
SmoLogFertRate = zeros(nt,n);

for k=2:58
    fspl = fit(logFerRate(2:14,1),logFerRate(2:14,k),'smoothingspline');
    SmoLogFertRate(:,(k-1)) = feval(fspl, u);

end
SmoLogFertRate

%*******************************************************************%
%*******************************************************************%

% Nessa etapa será feita uma análise de componentes principais funcional
% para as taxas em função das maturidades, onde serão usadas os funcionais
% discretizados (avaliados em pontos específicos), tomando como base o PCA
% descrito em Ramsay e Silverman (2005, p. 162 e 163). A base usada na
% expansão dos funcionais será a base de ondaletas.

% média cruzada (ao longo dos meses) das taxas
mu = mean(SmoLogFertRate,2)
% matriz com funcionais discretizados e centrados (com média zero)
X = (SmoLogFertRate - mu*ones(1,n))';

N = 3;  % nível usado na função wavedec
wname = 'db4';  % base do ondaletas usada nas decomposições
J = 84;
Xdec = zeros(n,J);

for ii = 1:n
    [Xdec(ii,:),Lw] = wavedec(X(ii,:),N,wname);
end

% versão discretizada do operador de variância-covariância
mV = (1/n)*(Xdec')*Xdec;

% Análise de autovalores e autovetores. As colunas da matriz B serão os
% coeficientes da decomposição das autofunções discretizadas, enquanto que
% os valores na diagonal de L serão seus autovalores correspondentes
[B,L] = eig(mV);

sprintf('%8.2f, ',diag(L(1:5,1:4)))  % maiores autovalores
% porcentagem da variabilidade explicada por cada autofunção
100*diag(L(1:4,1:4))/norm(diag(L),1)

% recuperando autofunções avaliadas nos pontos do vetor u
h1_hat =  waverec(B(:,1),Lw,wname);
h2_hat =  waverec(B(:,2),Lw,wname);
h3_hat =  waverec(B(:,3),Lw,wname);




vbeta1 = X*h1_hat;
vbeta2 = X*h2_hat;
vbeta3 = X*h3_hat;
%subplot(1,3,1); plot(2470:2500,vbeta1)
%subplot(1,3,2); plot(2470:2500,vbeta2)
%subplot(1,3,3); plot(2470:2500,vbeta3)

subplot(1,4,1); plot(u,mu);
subplot(1,4,2); plot(u,h1_hat);
subplot(1,4,3); plot(u,h2_hat);
subplot(1,4,4); plot(u,h3_hat);

%*******************************************************************%
%*******************************************************************%

%plot(diag(L(1:5,1:5)))

p = 5;  % lag máximo utilizado
N = 3;  % nível usado na função wavedec
wname = 'db4';  % base do ondaletas usada nas decomposições
% número de coefientes da base de ondaletas usada para o nt utilizado
J = 84;
% matriz de cujas linhas terão a decomposição em ondaletas
% de cada funcional observado (linha da matriz X)
Xdec = zeros(J,n);

% decomposição de cada funcional dos n anos
for ii = 1:n
    [Xdec(:,ii),Lw] = wavedec(SmoLogFertRate(:,ii),N,wname);
end

mu_dec = mean(Xdec,2);

C = Xdec - mu_dec*ones(1,n);

C1 = C(:,1:(n-p));
D1 = zeros(n-p,n-p);
for k=1:p
    D1 = D1 + C(:,(k+1):(n-p+k))'*C(:,(k+1):(n-p+k));
end
D = C1*D1*C1'/((n-p)^2);
[B,L] = eig(D);

diag(L(1:5,1:5))'
round(100*real(diag(L(1:5,1:5)))/norm(diag(L)),4)'

Nboot = 301;  % número de réplicas bootstrap
alpha = .05;  % nível de significância do teste bootstrap
% realizando os testes bootstrap para estimar a dimensão do processo
d0 = 1;
mPvalues = ones(1,8);  % valores-p dos testes dos 8 maiores autovalores

% Usando uma semente para o teste bootstrap
rng(2018);
while d0<=8
    % simularemos um processo de dimensão d0 e assim obtemos a distribuição
    % empírica bootstrap do (d0+1)-ésimo maior autovalor de D usando a
    % função DimEst_wavestrap
    d_boot = DimEst_wavestrap( Xdec, Nboot, B(:,[1:d0]), p);
    % valor-p para o (d0+1)-ésimo maior autovalor de D
    mPvalues(1,d0) = sum(d_boot>L(d0+1,d0+1))/Nboot;
    % verificamos se a hipótese nula de que esse autovalor é zero não é
    % rejeitada
    if (mPvalues(1,d0)>alpha)
        % se a hipótese não for rejeitada para d0+1, então a dimensão
        % selecionada será d0
        vDimSel = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end

mPvalues
vDimSel

mu1 =  waverec(mu_dec,Lw,wname);
% recuperando autofunções avaliadas nos pontos do vetor u
h1_hat =  waverec(B(:,1),Lw,wname);
h2_hat =  waverec(B(:,2),Lw,wname);
h3_hat =  waverec(B(:,3),Lw,wname);
h4_hat =  waverec(B(:,4),Lw,wname);

%*******************************************************************%
%*******************************************************************%

% Nessa simulação usamos o terceiro tipo de bootstrap, em que uma
% limiarização é aplicada na função média e nas autofunções para cirar um
% resíduo que será usado para formar as amostras bootstrap

Nboot = 301;  % número de réplicas bootstrap
alpha = .05;  % nível de significância do teste bootstrap
% realizando os testes bootstrap para estimar a dimensão do processo
d0 = 1;
mPvalues = ones(1,8);  % valores-p dos testes dos 8 maiores autovalores

% Usando uma semente para o teste bootstrap
while d0<=8
    % simularemos um processo de dimensão d0 e assim obtemos a distribuição
    % empírica bootstrap do (d0+1)-ésimo maior autovalor de D usando a
    % função DimEst_wavestrap
    d_boot = DimEst_wavestrap_thresh(Xdec, Nboot, B(:,[1:d0]), p, N, wname);
    % valor-p para o (d0+1)-ésimo maior autovalor de D
    mPvalues(1,d0) = sum(d_boot>L(d0+1,d0+1))/Nboot;
    % verificamos se a hipótese nula de que esse autovalor é zero não é
    % rejeitada
    if (mPvalues(1,d0)>alpha)
        % se a hipótese não for rejeitada para d0+1, então a dimensão
        % selecionada será d0
        vDimSel = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end

mPvalues
vDimSel

subplot(1,5,1); plot(u,mu);
ylabel('$\mu(x)$','Interpreter','latex');xlabel('Log tempo');
subplot(1,5,2); plot(u,mu1);
ylabel('$h_1(x)$','Interpreter','latex');xlabel('Log tempo');
subplot(1,5,3); plot(u,h2_hat);
ylabel('$h_2(x)$','Interpreter','latex');xlabel('Log tempo');
subplot(1,5,4); plot(u,h3_hat);
ylabel('$h_3(x)$','Interpreter','latex');xlabel('Log tempo')
subplot(1,5,5); plot(u,h4_hat);
ylabel('$h_4(x)$','Interpreter','latex');xlabel('Log tempo');

plot(SmoLogFertRate(2:14,:))

%*******************************************************************%
%*******************************************************************%

% Nessa simulação usamos o segundo tipo de bootstrap, em que uma
% limiarização é aplicada nos funcionais observados antes de realizar um
% procedimento de bootstrap de Bathia et al (2010)

p = 5;  % lag máximo utilizado
N = 3;  % nível usado na função wavedec
wname = 'db4';  % base do ondaletas usada nas decomposições
% número de coefientes da base de ondaletas usada para o nt utilizado
J = 84;
% matriz de cujas linhas terão a decomposição em ondaletas
% de cada funcional observado (linha da matriz X)
Xdec = zeros(J,n);

% decomposição de cada funcional dos n anos
for ii = 1:n
    % As linhas de Xdec vão armazenar as decomposições das linhas de X,
    % i.e., dos funcionais. Lw vai armazenar os númeroes de coeficientes em
    % cada nível das decomposições, sendo seu último elemente o número de
    % pontos das colunas
    [~,Xdec(:,ii),Lw] = wden(SmoLogFertRate(:,ii),'sqtwolog','h','mln',N,wname);
end

% esse vetor terá a média de um determinado coeficiente
% ao longo dos n dias
mu_dec = mean(Xdec,2);

% matriz com os desvios dos coeficientes em relação à
% média dos coeficientes num mesmo dia
C = Xdec - mu_dec*ones(1,n);

% agora será obtida a matriz formada pelos coeficientes das
% decomposições das observações dos funcionais. 
C1 = C(:,1:(n-p));
D1 = zeros(n-p,n-p);
for k=1:p
    D1 = D1 + C(:,(k+1):(n-p+k))'*C(:,(k+1):(n-p+k));
end
D = C1*D1*C1'/((n-p)^2);

% após obter D, calculamos os autovetores e autovalores da mesma.
% as colunas de B terão os autovetores de D, que são os coeficientes de
% ondaletas das autofunções associadas aos autovalores na diagonal de L
[B,L] = eig(D);

Nboot = 301;  % número de réplicas bootstrap
alpha = .05;  % nível de significância do teste bootstrap
% realizando os testes bootstrap para estimar a dimensão do processo
d0 = 1;
mPvalues = ones(1,8);  % valores-p dos testes dos 8 maiores autovalores

while d0<=8
    % simularemos um processo de dimensão d0 e assim obtemos a distribuição
    % empírica bootstrap do (d0+1)-ésimo maior autovalor de D usando a
    % função DimEst_wavestrap
    d_boot = DimEst_wavestrap( Xdec, Nboot, B(:,[1:d0]), p);
    % valor-p para o (d0+1)-ésimo maior autovalor de D
    mPvalues(1,d0) = sum(d_boot>L(d0+1,d0+1))/Nboot;
    % verificamos se a hipótese nula de que esse autovalor é zero não é
    % rejeitada
    if (mPvalues(1,d0)>alpha)
        % se a hipótese não for rejeitada para d0+1, então a dimensão
        % selecionada será d0
        vDimSel = d0;
        d0 = 9;
    end
    d0 = d0 + 1;
end

mPvalues
vDimSel

