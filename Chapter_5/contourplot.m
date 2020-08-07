%Second Generation Wavelets. Mapping the contour plot. 
%% This chunk sets up u
clear
tend=5.2;
g=10;
n=2^9; %grid points
b=2*pi; %length of x axis
delx= b/n; %width of space step
delt=0.1*delx;
w=length(0:delt:tend);
visc=delx^2/8
x= 0:delx:b-delx; %adds delx each time and specifies grid points
uinit=zeros(1,n); %preallocating u
for i=1:n
    uinit(i)= sin(x(i));
end
%[t,u]=rk4(@ode45try2,0, tend, uinit, w); 
u(1,:)=uinit;
len=length(x);
for p=1:w
    [t,thing]=rk4try2(@rk4setup,delt*(p-1), delt*p, u(p,:), 1,len,1);
    u(p+1,:)=thing(2,:);
end

%% This chunk runs the contour plot stuff
for n=w
s = [u(n,:)]; 
len = length(s); %number of grid points

lev   = 4;
nbcol = 100;
App=zeros(lev, (length(s))/2); 
Dt=zeros(lev,[length(s)]/2); 
cfd=zeros(lev,[length(s)]);

%perform decomposition
[App(1,:),Dt(1,:)]=waveinter(s,1,0);

for i=2:lev
     Ex = App(i-1,1:((length(x))/(2^(i-1))));
    [App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = waveinter(Ex, 1,0);
end
for i=1:lev
    d=Dt(i,1:(len/(2^i)));
    d=d(:)';
    d=d(ones(1,2^i),:);
    cfd(i,:)=wkeep1(d(:)',len);
end

cfd =  cfd(:);
% I = find(abs(cfd)<sqrt(eps));
% cfd(I) = zeros(size(I));
cfd    = reshape(cfd,lev,len);
cfd = wcodemat(cfd,nbcol,'row'); %put in ,'row' after nbcol if you want to scale by row


h211 = subplot(211);
h211.XTick = [];
plot(s,'r');grid on;
xlim([0,len])
set(gca, 'XTick', [0:0.1:1]*len, 'XTickLabel', [0:0.1:1]*2)

title('Burgers Equation');
ax = gca;
%ax.XLim = [1 length(s)];
%ax.XLim = [0 2]
ax.YLim = [-1 1]
% % 
% indices = find(abs(cfd)<80); %delete values below a certain magnitude
% cfd(indices) = 0;  

subplot(212);
%colormap(flipud(gray));
colormap(bone(80))

image(cfd);
tics = 1:lev;
labs = int2str((fliplr(tics))');
ax = gca;
ax.YTickLabelMode = 'manual';
ax.YDir = 'normal';
ax.Box = 'On';
ax.YTick = tics;
ax.YTickLabel = labs;

set(gca, 'XTick', [0:0.1:1]*len, 'XTickLabel', [0:0.1:1]*2)
%set(gca, 'XTick', [0:0.1:20]*200, 'XTickLabel', [0:0.1:20]*2) %makes axis of contour as x/2pi

title('Wavelet coefficient magnitude');
ylabel('Resolution Level');

p = colorbar('southoutside')
set(p,'ylim',[1 100])


mov(n)=getframe(figure(1));

end

% v = VideoWriter('contourpresentation.avi');
% v.FrameRate = 110;  % Default 30
% v.Quality = 100;    % Default 75
% open(v)
% writeVideo(v,mov)
% close(v)    