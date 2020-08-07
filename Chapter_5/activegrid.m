%% this sets approximation to zero. so this is just for visualization purposes of the adaptive grid.
% Look at the activegridcalc.m file for the one that just sets the wavelets
% to zero, but keeps track of which ones should be zeroed out.
function [App, Dt]=activegrid(App, Dt, eps, lev)
%This function spits out the Dt, App, and y1 (finest resolution) structures
%for the active grid of a function with certain threshold on the wavelets.

I2=find(abs(Dt')>eps); %finds points on the grid which are currently above the threshold. counts down column usually, but we want it to count down row.

len=length(App(1,:))*2;

%Security/Safety zone

I3=I2;%store values in I3
for k=1:len*lev
if ismember(k,I2)==1 %if we are in the I2 matrix
    for i=1:lev
        if k>(len/2)*(i-1) && k<=(len/2)*(i) %determine which line the point is on.
            row=i;
        end
    end
        m=k-(len/2)*(row-1); %gives position on line
        if row ~=1 %if we aren't on the first row
            I3(length(I3)+1)=k-m-(len/2)+ 2*m; 
            I3(length(I3)+1)=k-m-(len/2)+2*m -1;
        end
end
end
%Remember, we have to transpose it because we counted the transposed way.    
B=Dt';
I4=0;
for v=1:lev*(len/2)
    if ismember(v,I3)==0
        I4(length(I4)+1)=v; %store all wavelet values not in safety zone
    end
end
I4=I4(2:end)';
B(I4)=zeros(size(I4));
Dt=B';


%Modify the app matrix to find those at coarsest level and in the perfect
%reconstruction zone
%find wavelets at coarsest level
for i=1:(len/(2^lev))
I5(i)=(lev-1)*(len/2)+i;
end

%find those in perfect reconstruction zone above each active wavelet in I2.
%also add in grid pts directly above active wavelets (line 61 and 66). idk
%if should??????
for k=1:len*lev
if ismember(k,I3)==1 %if we are in the I2 matrix
    for i=1:lev
        if k>(len/2)*(i-1) && k<=(len/2)*(i) %determine which line the point is on.
            row=i;
        end
    end
        m=k-(len/2)*(row-1); %gives position on line
        if row~=1 && m~=len/2^row %if not the last wavelet on a line and not on the first line
           I5(length(I5)+1)=k-m-(len/2-1) +2*m-2; %gets the scaling function to the left of active wavelet
            I5(length(I5)+1)=k-m-(len/2-1) +2*m-1; %gets the scaling function directly above the active wavelet ???????? SHOULD I
           I5(length(I5)+1)=k-m-(len/2-1)+2*m; %line above and to the right
        end
        if row~=1 && m==len/2^row %if not on the first line, but are the last wavelet on a line
           I5(length(I5)+1)=k-m-(len/2 -1) + 2*m -2;
            I5(length(I5)+1)=k-m-(len/2-1) +2*m-1; %gets the scaling function directly above the active wavelet ?????? SHOULD I
           I5(length(I5)+1)=k-m-(len/2 -1); %gets first scaling function point on line above
        end  
    end
end


%Remember, we have to transpose it because we counted the transposed way. 
%Make everything in the App matrix 0 unless it is in the perfect
%reconstruction zone or the coarsest level
C=App';
I7=0;
for v=1:lev*(len/2)
    if ismember(v,I5)==0
        I7(length(I7)+1)=v; %store all scaling function values not in perfect reconstruction zone or in coarsest level
    end
end
I7=I7(2:end)';
C(I7)=zeros(size(I7));
App=C';
end
