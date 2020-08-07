%this spits out the Dt and App matrices and a tracker matrix which shows us which points to
%keep (1=keep, 0 = don't keep).
function [App, Dt,tr]=activegridcalc(App, Dt, eps, lev)
[A, D]=activegrid(App, Dt, eps, lev);
len=length(App(1,:))*2;

tr=zeros(lev,len);
Dt1=zeros(lev,len);
App1=zeros(lev,len);
for i=1:lev
    x=D(i,1:len/(2^i));
    Dt1(i,1+(2^(i-1)):(2^i):end)=x;
    x1=A(i,1:len/(2^i));
    App1(i,1:(2^i):end)=x1;
    tr(i,:)=Dt1(i,:)+App1(i,:);
end
I=find(abs(tr)>0);
I2=find(abs(tr)==0);
%tr(I2)=NaN(size(I2));
tr(I2)=zeros(size(I2));
tr(I)=ones(size(I));
App=A;
Dt=D;
end