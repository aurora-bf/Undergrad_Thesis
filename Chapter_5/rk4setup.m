function udot=rk4setup(t,u,len,m)
%len is the original signal length
%m is the row we are on.
nn=len/(2^(m-1));
b=2*pi; %length of x axis
delx= b/nn; %width of space step
%visc=delx^1.2; %viscosity coefficient
visc=(0.0123^2)/8;
for i=1:nn
    if i==1
%         udot(i,1)= -u(1) * (u(2)-u(n))/(2*delx) + visc * (u(2)-2*u(1)+u(n))/delx^2;
        udot(i,1)= -u(1) * (u(2)-u(nn-1))/(2*delx) + visc * (u(2)-2*u(1)+u(nn-1))/delx^2;
    elseif i<nn
        udot(i,1) = -u(i)*(u(i+1)-u(i-1))/(2*delx) + visc*(u(i+1)-2*u(i)+u(i-1))/delx^2;
    elseif i==nn
          udot(i,1)= -u(nn)*(u(1)-u(nn-1))/(2*delx) + visc*(u(1)-2*u(nn)+u(nn-1))/delx^2;
%           udot(i,1)=u(1,1)
    end
end
end

        

%     udot(1) =  u(1) * (u(2)-u(n))/(2*delx) + visc * (u(2)-2*u(1)+u(n))/delx^2;
%     
%     for i=2:n-1
%         udot(i) = u(i)*(u(i+1)-u(i-1))/(2*delx) + visc*(u(i+1)-2*u(i)+u(i-1))/delx^2 ; %from ui=ui+udot(ui)
%     end
%    
%     udot(n)= (u(n)*(u(1)-u(n-1)))/(2*delx) + visc*(u(1)-2*u(n)+u(n-1))/delx^2 ;



% put ode45 in a separate file
% [t,u]=ode45(@ex6_ode,0:delt:tend,u,[],visc);