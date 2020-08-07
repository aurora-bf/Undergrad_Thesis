
function [A,D]= waveinter(y,m,e)
if mod(length(y),2)~=0
    msg = 'Error occurred. Length of initial function must be even.';
    error(msg)
end
%A is the approximation, D is the detail, a3 is the coefficient of lifting
%to preserve mean
%m=order of interpolation
if m==1 %linear interpolation
d=zeros(1,length(y));
a=zeros(1,length(y));
%forward transform
for i=1:length(y)
if (mod(i,2)==0) && (i~=length(y));%if even and not an end point
    d(i)=0.5*(y(i)-0.5*(y(i-1)+y(i+1))); %define wavelet definition
end
if mod(i,2)==0 && i==length(y); %if even and is an end point
    d(i)=0.5*(y(i)-0.5*(y(i-1)+y(1)));
end
end 
I=find(abs(d)<e);
d(I)=0;
%define approximation
for i=1:length(y)
if mod(i,2)~=0 && i~=1; %if odd and not an end point
    a(i)=y(i)+0.5*(d(i-1)+d(i+1));
end
if mod(i,2)~=0 && i==1
    a(i)=y(i)+0.5*(d(i+1)+d(length(y)));
end
end
A=a(1:2:end);
D=d(2:2:end);
end

%quadratic interpolation
if m==2
d=zeros(1,length(y));
a=zeros(1,length(y));
z=zeros(1,length(y));
%define the z values (where z(i) is the deriv of the interpolant at pt i)
for i=1:length(y)
   if (mod(i,2)~=0) && (i~=length(y)-1) && (i~=1);%if even and not an end point
     z(i)=y(i+2)-y(i-2); %define wavelet definition
end
if mod(i,2)~=0 && i==length(y)-1; %if even and is an end point
    z(i)=y(1)-y(i-2);
end
if mod(i,2)~=0 && i==1; %if even and is an end point
    z(i)=y(i+2)-y(length(y)-1);

end
end
%define wavelets
for i=1:length(y)
if (mod(i,2)==0) && (i~=length(y));%if even and not an end point
d(i)=0.5*(y(i)-(0.75*z(i-1)+(1/4)*z(i+1)+y(i-1)));
end
if mod(i,2)==0 && i==length(y); %if even and is an end point
d(i)=0.5*(y(i)-(0.75*z(i-1)+(1/4)*z(1)+y(i-1)));
end
end
I=find(abs(d)<e);
d(I)=0;
%define approximation
for i=1:length(y)
if mod(i,2)~=0 && i~=1; %if odd and not an end point
    a(i)=y(i)+0.5*(d(i-1)+d(i+1));
end
if mod(i,2)~=0 && i==1
    a(i)=y(i)+0.5*(d(i+1)+d(length(y)));
end
end
A=a(1:2:end);
D=d(2:2:end);
end

       

    
    
end
    
    
    
    
    
    
    


%quadratic interpolation LAGRANGE :(
% if m==2
%     d=zeros(1,length(y));
% a=zeros(1,length(y));
% %forward transform
% for i=1:length(y)
% if (mod(i,2)==0) && (i~=length(y)) && (i~=length(y)-2);%if even and not an end point or a pt 3 from the end
%     d(i)=0.5*(y(i)-(3/8)*y(i-1)-(3/8)*y(i+1)-(1/8)*y(i+3)); %define wavelet definition
% end
% if mod(i,2)==0 && i==length(y); %if even and is an end point
%     d(i)=0.5*(y(i)-(3/8)*y(i-1)-(3/8)*y(1)-(1/8)*y(3));
% end
% if mod(i,2)==0 && i==length(y)-2
%     d(i)= 0.5*(y(i)-(3/8)*y(i-1)-(3/8)*y(i+1)-(1/8)*y(1))
% end 
% end
% 
% %define the coefficient a
% a1=0
% a2=0
% for j=1:length(y)
%   a1=4*y(j)*(-1)^(j-1) +a1;
%   if mod(j,2)~=0
%       a2=7*y(j) +a2;
%   end
%   if mod(j,2)==0
%       a2=-8*y(j) +a2;
%   end
% end
% a3=a1/a2
% 
% for i=1:length(y)
% if mod(i,2)~=0 && i~=1; %if odd and not an end point
%     a(i)=y(i)+a3*(d(i-1)+d(i+1));
% end
% if mod(i,2)~=0 && i==1
%     a(i)=y(i)+a3*(d(i+1)+d(length(y)));
% end
% end
% A=a(1:2:end);
% D=d(2:2:end);
% end











