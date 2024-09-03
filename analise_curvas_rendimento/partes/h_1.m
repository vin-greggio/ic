logFerRate = csvread('parte1_f.csv');
logFerRate2 = csvread('parte2_f.csv');
logFerRate3 = csvread('parte3_f.csv');
logFerRate4 = csvread('parte4_f.csv');
logFerRate5 = csvread('parte5_f.csv');
logFerRate6 = csvread('parte6_.csv');
logFerRate7 = csvread('parte7_.csv');
logFerRate8 = csvread('parte8_.csv');
logFerRate9 = csvread('parte9_.csv');
logFerRate10 = csvread('parte10_.csv');
logFerRate11 = csvread('parte11_.csv');
logFerRate12 = csvread('parte12_.csv');
logFerRate13 = csvread('parte13_.csv');
logFerRate14 = csvread('parte14_.csv');
logFerRate15 = csvread('parte15_.csv');
logFerRate16 = csvread('parte16_.csv');

n = 54;
nt = 64;
u = linspace(0,5.886104031450156,nt);

SmoLogFertRate1 = zeros(nt,n);

for k=2:55
    fspl = fit(logFerRate(2:14,1),logFerRate(2:14,k),'smoothingspline');
    SmoLogFertRate1(:,(k-1)) = feval(fspl, u);
end
n = 60;

SmoLogFertRate2 = zeros(nt,n);
SmoLogFertRate3 = zeros(nt,n);
SmoLogFertRate4 = zeros(nt,n);
SmoLogFertRate5 = zeros(nt,n);
SmoLogFertRate6 = zeros(nt,n);
SmoLogFertRate7 = zeros(nt,n);
SmoLogFertRate8 = zeros(nt,n);
SmoLogFertRate9 = zeros(nt,n);
SmoLogFertRate10 = zeros(nt,n);
SmoLogFertRate11 = zeros(nt,n);
SmoLogFertRate12 = zeros(nt,n);
SmoLogFertRate13 = zeros(nt,n);
SmoLogFertRate14 = zeros(nt,n);
SmoLogFertRate15 = zeros(nt,n);
SmoLogFertRate16 = zeros(nt,n);


for k=2:61
    fspl = fit(logFerRate2(2:14,1),logFerRate2(2:14,k),'smoothingspline');
    SmoLogFertRate2(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate3(2:14,1),logFerRate3(2:14,k),'smoothingspline');
    SmoLogFertRate3(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate4(2:14,1),logFerRate4(2:14,k),'smoothingspline');
    SmoLogFertRate4(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate5(2:14,1),logFerRate5(2:14,k),'smoothingspline');
    SmoLogFertRate5(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate6(2:14,1),logFerRate6(2:14,k),'smoothingspline');
    SmoLogFertRate6(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate7(2:14,1),logFerRate7(2:14,k),'smoothingspline');
    SmoLogFertRate7(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate8(2:14,1),logFerRate8(2:14,k),'smoothingspline');
    SmoLogFertRate8(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate9(2:14,1),logFerRate9(2:14,k),'smoothingspline');
    SmoLogFertRate9(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate10(2:14,1),logFerRate10(2:14,k),'smoothingspline');
    SmoLogFertRate10(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate11(2:14,1),logFerRate11(2:14,k),'smoothingspline');
    SmoLogFertRate11(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate12(2:14,1),logFerRate12(2:14,k),'smoothingspline');
    SmoLogFertRate12(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate13(2:14,1),logFerRate13(2:14,k),'smoothingspline');
    SmoLogFertRate13(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate14(2:14,1),logFerRate14(2:14,k),'smoothingspline');
    SmoLogFertRate14(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate15(2:14,1),logFerRate15(2:14,k),'smoothingspline');
    SmoLogFertRate15(:,(k-1)) = feval(fspl, u);

end
for k=2:61
    fspl = fit(logFerRate16(2:14,1),logFerRate16(2:14,k),'smoothingspline');
    SmoLogFertRate16(:,(k-1)) = feval(fspl, u);

end

valores_lista = []
h_lista = {}
p_lista = []
smo_lista = {SmoLogFertRate1;SmoLogFertRate2;SmoLogFertRate3;SmoLogFertRate4;SmoLogFertRate5;SmoLogFertRate6;SmoLogFertRate7;SmoLogFertRate8;SmoLogFertRate9;SmoLogFertRate10;SmoLogFertRate11;SmoLogFertRate12;SmoLogFertRate13;SmoLogFertRate14;SmoLogFertRate15;SmoLogFertRate16}

for x=1:16
    SmoLogFertRate = smo_lista{x,1};
    if x==1
        n=54;
    else
        n=60;
    end
    mu = mean(SmoLogFertRate,2);
    X = (SmoLogFertRate - mu*ones(1,n))';
    
    p = 5;
    N = 3; 
    wname = 'db5'; 

    J = 88;
    Xdec = zeros(J,n);

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
    100*diag(L(1:4,1:4))/norm(diag(L),1)
    %diag(L(1:5,1:5))'
    mu =  waverec(mu_dec,Lw,wname);
    h1_hat =  waverec(B(:,1),Lw,wname);
    h2_hat =  waverec(B(:,2),Lw,wname);
    h3_hat =  waverec(B(:,3),Lw,wname);
    h4_hat =  waverec(B(:,4),Lw,wname);


    m_lista = [h1_hat,h2_hat,h3_hat,h4_hat];
    h_lista{end+1} = m_lista;

end



h_lista{1,16}(:,1) = -h_lista{1,16}(:,1);
h_lista{1,15}(:,1) = -h_lista{1,15}(:,1);
h_lista{1,14}(:,1) = -h_lista{1,14}(:,1);
h_lista{1,8}(:,1) = -h_lista{1,8}(:,1);
h_lista{1,2}(:,1) = -h_lista{1,2}(:,1);



listamu = [h_lista{1,1}(:,1),h_lista{1,2}(:,1),h_lista{1,3}(:,1),h_lista{1,4}(:,1),h_lista{1,5}(:,1),h_lista{1,6}(:,1),h_lista{1,7}(:,1),h_lista{1,8}(:,1),h_lista{1,9}(:,1),h_lista{1,10}(:,1),h_lista{1,11}(:,1),h_lista{1,12}(:,1),h_lista{1,13}(:,1),h_lista{1,14}(:,1),h_lista{1,15}(:,1),h_lista{1,16}(:,1)];

mesh([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],u,listamu)
