function [yt]=waveinterinv(A,D,m)
%1st order interpolation
if m==1
a = reshape([A; zeros(size(A))],[],1)'; %put the zeroes back into a
d = reshape([zeros(size(D)); D],[],1)'; %put the zeroes back into d
%Inverse transform
for i=1:length(a)
if mod(i,2)~=0 && i~=1; %if odd and not an end point
    yt(i)=a(i)-0.5*(d(i-1)+d(i+1));
end
if mod(i,2)~=0 && i==1
    yt(i)=a(i)-0.5*(d(i+1)+d(length(a)));
end
end

for i=1:length(a)
if (mod(i,2)==0) && (i~=length(a));%if even and not an end point
    yt(i)= 2*d(i)+0.5*(yt(i-1)+yt(i+1)); %reconstruction
end
if mod(i,2)==0 && i==length(a); %if even and is an end point
    yt(i)= 2*d(i)+0.5*(yt(i-1)+yt(1));
end
end
end


if m==2
a = reshape([A; zeros(size(A))],[],1)'; %put the zeroes back into a
d = reshape([zeros(size(D)); D],[],1)'; %put the zeroes back into d
%Inverse transform
for i=1:length(a)
if mod(i,2)~=0 && i~=1; %if odd and not an end point
    yt(i)=a(i)-0.5*(d(i-1)+d(i+1));
end
if mod(i,2)~=0 && i==1
    yt(i)=a(i)-0.5*(d(i+1)+d(length(a)));
end
end


%define new z functions
for i=1:length(a)
   if (mod(i,2)~=0) && (i~=length(a)-1) && (i~=1);%if even and not an end point
    zt(i)=yt(i+2)-yt(i-2); %define wavelet definition
end
if mod(i,2)~=0 && i==length(a)-1; %if even and is an end point
    zt(i)=yt(1)-yt(i-2);
end
if mod(i,2)~=0 && i==1; %if even and is an end point
    zt(i)=yt(i+2)-yt(length(a)-1);
end
end
 %recover the even functions
for i=1:length(a)
if (mod(i,2)==0) && (i~=length(a));%if even and not an end point
  yt(i)= 2*d(i)+0.75*zt(i-1)+(1/4)*zt(i+1)+yt(i-1); %reconstruction
end
if mod(i,2)==0 && i==length(a); %if even and is an end point
yt(i)= 2*d(i)+0.75*zt(i-1)+(1/4)*zt(1)+yt(i-1);
end
end
end




%2nd order interpolation FAILED LAGRANGE
% if m==2
% a = reshape([A; zeros(size(A))],[],1)'; %put the zeroes back into a
% d = reshape([zeros(size(D)); D],[],1)'; %put the zeroes back into d
% %Inverse transform
% for i=1:length(a)
% if mod(i,2)~=0 && i~=1; %if odd and not an end point
%     yt(i)=a(i)-a3*(d(i-1)+d(i+1));
% end
% if mod(i,2)~=0 && i==1
%     yt(i)=a(i)-a3*(d(i+1)+d(length(a)));
% end
% end
% 
% for i=1:length(a)
% if (mod(i,2)==0) && (i~=length(a)) && (i~=length(a)-2);%if even and not an end point or a pt 3 from the end
%     yt(i)=2*d(i)+(3/8)*yt(i-1)+(3/8)*yt(i+1)+(1/8)*yt(i+3); %define wavelet definition
% end
% if mod(i,2)==0 && i==length(a); %if even and is an end point
%     yt(i)=2*d(i)+(3/8)*yt(i-1)+(3/8)*yt(1)+(1/8)*yt(3);
% end
% if mod(i,2)==0 && i==length(a)-2
%     yt(i)= 2*d(i)+(3/8)*yt(i-1)+(3/8)*yt(i+1)+(1/8)*yt(1)
% end 
% end
% end
end