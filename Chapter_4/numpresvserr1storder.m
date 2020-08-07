for b=1:2
j=16;
n=14;%number of levels of resolution analyzed. 
xend=4;
xbeg=-4;
lenneeded =length(xbeg:(xend-xbeg +(1/2^(j)))/(2^(j)):xend); 
x=linspace(xbeg,xend,lenneeded);
nbcol=100;
x1=linspace(xbeg,xend,lenneeded/2);

% if b==1
% y=x.^2
% end

if b==1
y1=tanh(x1);
y2=-tanh(x1);
y=[y1 y2];
end

if b==2
y=gaussmf(x,[0.5 0]);
end





p=66; %number of epsilons tested
%p=60
q=1 %order of interpolation tested


for m=1:q
for k=1:p %choosing how many different epsilons to test
   th=((k/10)-8.1)
    thres(k,:)=th;
e=10^th; %setting an epsilon
% e=k*10^-4;
ep(k,:)=e; %vary threshold and keep track

len=length(x); 
App=zeros(n, (length(x))/2); 
Dt=zeros(n,[length(x)]/2); 

[App(1,:),Dt(1,:)]=waveinter(y,m,e); %first decomposition

for i=2:n
    Ex = App(i-1,1:((length(x))/(2^(i-1))));
[App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = waveinter(Ex, m,e);
end
I2=find(abs(Dt)>0);
numpres(k,:) = prod(size(I2));



yt= zeros(n, (length(x)));
yt(1, 1:((len/(2^(n-1)))))=waveinterinv(App(n,1:((len/(2^n)))),Dt(n,1:(len/(2^n))),m); %start with y1, reconstruct

for i=2:n
    yt(i, 1:(len/(2^(n-i))))=waveinterinv(yt(i-1, 1:(len/(2^((n-i+1))))),Dt((n+1-i),1:(len/(2^(n-i+1)))),m); %reconstruct up all levels
end

% err(k,:) = max(abs(yt(n,:)-y)); %check perfect reconstruction
err(k,:)=norm(yt(n,:)-y,2)
end
%set zero points to one so we can see them
% I3=find(numpres==0)
% numpres(I3)=ones(size(I3))
end

figure(1)
loglog(numpres,err,'-')
title('Number of Preserved Collocation Points versus Error')
ylabel('L^2 Norm of Error')
xlabel('Number of Preserved Collocation Points')
hold on;
% legendInfo{m} = [;'Polynomial Interpolation of Order' num2str(m)]
% legend(legendInfo)
% legend('Location','southwest')
grad=polyfit(log(numpres),log(err),1)
slope(b)=grad(:,1)
 
end
loglog(numpres,numpres.^(-2),'--k')
 legend('tanh(x)','gaussmf(x)','Reference line with slope -2')
