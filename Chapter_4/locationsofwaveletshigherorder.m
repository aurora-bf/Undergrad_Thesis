clear
n=5; %number of levels of resolution analyzed. 
j=8; %keep j>n. determines length signal
%Max resolution level that this code works for is for where
%(length(x))/(2^n) is still an integer (for even lengtheds things)
xend=4;
xbeg=-4;
lenneeded =length(xbeg:(xend-xbeg +(1/2^(j)))/(2^(j)):xend); %having it go up by (xend-xbeg/(2^n)) ensures that length(x)/(2^n) is an integer
x=linspace(xbeg,xend,lenneeded);
% 
x1=linspace(xbeg,xend,lenneeded/2)
y1=tanh(x1);
y2=-tanh(x1);
y=[y1 y2];

nbcol=100;
%y= gaussmf(x,[0.5 0]);
%y=zeros(1,length(x));

%%
% m=3;
% len=length(x); %The length of the approx matrix is = to length(x). length of detail matrix is = to length(x)-1
% App=zeros(n, (length(x))/2); %initializes matrix for approx. Sets it to length of first level which is half of original
% Dt=zeros(n,[length(x)]/2); %initializes matrix for detail
% LS=liftwave('lazy');
% ElimLiftStep = {'d',[-1/m],0}; 
% LSNalmost=addlift(LS,ElimLiftStep,'end');
% elsprimal = {'p',[1/(2*m)], 0}; 
% LSN = addlift(LSNalmost,elsprimal,'end');
% %%
len=length(x); 
App=zeros(n, (length(x))/2); 
Dt=zeros(n,[length(x)]/2); 
e=0.0001
m=1

[App(1,:),Dt(1,:)]=waveinter(y,m,0); %first decomposition

for i=2:n
    Ex = App(i-1,1:((length(x))/(2^(i-1))));
[App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = waveinter(Ex, m,0);
end

I=find(abs(Dt)<e);
Dt(I)=NaN(size(I)); %set things below threshold to NaN so they don't plot

I2=find(abs(Dt)>=e);
Dt(I2)=zeros(size(I2)); %set things above the threshold to 0 so the points plot along the x axis


% %%
% 
% [App(1,:),Dt(1,:)]=lwt(y,LSN); %first decomposition
% 
% 
% %Below is the code for subsequent decompositions. Note that the size of
% %each successive level is the previous level/2
% for i=2:n
%     Ex = App(i-1,1:((length(x))/(2^(i-1))));
% [App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = lwt(Ex, LSN);
% end
% e=0.01;
% I=find(abs(Dt)<e);
% Dt(I)=NaN(size(I)); %set things below threshold to NaN so they don't plot
% 
% I2=find(abs(Dt)>=e);
% Dt(I2)=zeros(size(I2)); %set things above the threshold to 0 so the points plot along the x axis
% 
% 
% % for i=1:n
% % %xodd= x(2:2^i:end)
% % xodd=linspace(xbeg,xend,len/(2^i))
% % plot(xodd, (n-(i-1))/(n-1) -(1/(n-1)) + Dt(i,1:(len/2^i)),'.r','MarkerSize', 10) %plots show ehere wavelets are preserved at each resolution number
% % hold on;
% % plot(x,y)
% % end
%%
% xbeg=-4-2^(-n)%need to start so the first wavelet is shifted off by one step at the highest scale of resolution
store=zeros(n,len/2)


for i=1:n
yyaxis right
xodd=linspace(xbeg+((xend-xbeg)/[2^(j)-1])*(2^(i-1)),xend-((xend-xbeg)/[2^(j)-1])*(2^(i-1)-1),len/(2^i));%((xend-xbeg)/[2^(j)-1]) is the step size on the finest grid. (see 22 march page 2)
plot(xodd, (n-i) + 1 + Dt(i,1:((len/2^i))),'x','MarkerSize', 10) %plots show where wavelets are preserved at each resolution number
hold on;
ylabel('Resolution Level')
yticks([linspace(1,n,n)]);
store(i,1:(len/2^i))=xodd; %this is just to check and make sure everything lines up
end
yyaxis left
plot(x,y)
ylabel('Function value')
hold on; 





