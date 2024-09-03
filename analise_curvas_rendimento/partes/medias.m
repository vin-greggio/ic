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

SmoLogFertRate = zeros(nt,n);

for k=2:55
    fspl = fit(logFerRate(2:14,1),logFerRate(2:14,k),'smoothingspline');
    SmoLogFertRate(:,(k-1)) = feval(fspl, u);
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


mu = mean(SmoLogFertRate,2)
mu2 = mean(SmoLogFertRate2,2)
mu3 = mean(SmoLogFertRate3,2)
mu4 = mean(SmoLogFertRate4,2)
mu5 = mean(SmoLogFertRate5,2)
mu6 = mean(SmoLogFertRate6,2)
mu7 = mean(SmoLogFertRate7,2)
mu8 = mean(SmoLogFertRate8,2)
mu9 = mean(SmoLogFertRate9,2)
mu10 = mean(SmoLogFertRate10,2)
mu11 = mean(SmoLogFertRate11,2)
mu12 = mean(SmoLogFertRate12,2)
mu13 = mean(SmoLogFertRate13,2)
mu14 = mean(SmoLogFertRate14,2)
mu15 = mean(SmoLogFertRate15,2)
mu16 = mean(SmoLogFertRate16,2)

listamu = [mu,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10,mu11,mu12,mu13,mu14,mu15,mu16]

mesh([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],exp(u),listamu)

