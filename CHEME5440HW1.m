%Sam Furness
%2/3/19
%Written for HW1 CHEME5440

%Define Variables:
kEj=0.0137;%s^-1
RXT=19.8;%uM
Gj=6.196;%uM
kmin=0.1;%s^-1
kplus=5;%uM^-1s^-1
kI=1/42;%s^-1
kA=0;%s^-1
KXj= (kmin+kI)/kplus;
tau= (kA+kEj)/kI;
W1=0.26;
W2=300;
n=1.5;
K=0.3;
kd=log(2)/120;%s^-1
mu=1.14/3600;%s^-1

I=linspace(0.0001,10,100000);%mM
for i=1:length(I)
    fI(i)=(I(i)^n)/(K^n+I(i)^n);
end

for i=1:length(I)
    mj(i)=(((kEj*RXT*Gj)/(KXj*tau+Gj*(tau+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
end
    
logI=log10(I);
logmj=log10(mj);
semilogx(I,mj)
%plot(logI,mj);
%plot(I,logmj);
title('mRNA versus Inducer Concentration')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
    
    
    
    
