function y = IM_candidate(eta,beta,IM_end)
%{
This is the code for the identification of IM level candidates for the 
efficient fragility conversion method by Pang and Wang (202X):

Reference: 
Pang Y., and Wang X. (202X) “Cloud-IDA-MSA conversion of fragility curves
for efficient and high-fidelity resilience assessment? Journal of
Structural Engineering (ASCE), under review (STENG-9556R1)

Input:
   eta - vector of fragility median for multiple damage states (e.g. 4)
  beta - vector of fragility dispersion for multiple damage states (e.g. 4)
IM_end - scalar of the upper boundary of the considered IM range

Output:
     y - identified IM level candidates (vector)
%}

%% Define Cut-Off Probabilities for the Lower and Upper IM Boxes (Avoid Zero or Unity(1) Probabilities)
% Recommended values: 0.01 and 0.99 for the lower and upper ones,respectively.
% Users can modify to achieve more or less IM level candidates.
pf1=0.01; pf2=0.99;

%% Develop Original Fragility Curves and Determine the IM Boxes
im=linspace(0.001,IM_end,500);
hold on
for i=1:4
    pfc(:,i)= normcdf((log(im/eta(i)))/beta(i));
    im1(1,i)  =exp(norminv(pf1,0,1)*beta(i))*eta(i);
    im1(2,i+1)=exp(norminv(0.5,0,1)*beta(i))*eta(i);
    im2(1,i)  =exp(norminv(0.5,0,1)*beta(i))*eta(i);
    im2(2,i+1)=exp(norminv(pf2,0,1)*beta(i))*eta(i);
    plot(im,pfc(:,i),'k','linewidth',1);
end

%% Identify (Judge) Overlaps of the IM Boxes for Different Damage States (DS)
% DS1 & DS2: overlap?
DS=1;
if im1(1,DS+1)<=im1(2,DS+1)
    IM_can(DS,1)=0.5*(im1(1,DS+1)+im1(2,DS+1));
else
    IM_can(DS,1)=0.5*(im1(1,DS+1)+im1(1,DS));
end

if im2(1,DS+1)<=im2(2,DS+1)
    IM_can(DS,2)=0.5*(im2(1,DS+1)+im2(2,DS+1));
end

% DS2 & DS3: overlap
DS=2;
if im1(1,DS+1)<=im1(2,DS+1)
    IM_can(DS,1)=0.5*(im1(1,DS+1)+im1(2,DS+1));
end
if im2(1,DS+1)<=im2(2,DS+1)
    IM_can(DS,2)=0.5*(im2(1,DS+1)+im2(2,DS+1));
end

% DS3 & DS4: overlap
DS=3.0;
if im1(1,DS+1)<=im1(2,DS+1)
    IM_can(DS,1)=0.5*(im1(1,DS+1)+im1(2,DS+1));
end
if im2(1,DS+1)<=im2(2,DS+1)
    IM_can(DS,2)=0.5*(im2(1,DS+1)+im2(2,DS+1));
end

% DS4
[n,k]=size(IM_can);
m=1;
for i=1:n
    for j=1:k
        if IM_can(i,j)~=0
            IM_c(m,1)=IM_can(i,j);
            m=m+1;
        end
    end
end
IM_c=sort(IM_c);
n=find(IM_c>=IM_end);
if min(size(n))~=0
    IM_c(n)=IM_end;
end
IM_c=unique(IM_c);
if IM_c(end,1)<IM_end
    IM_c(end+1,1)=IM_end;
end

%% Draw the IM Level Candidates
IM_subscript={'a';'b';'c';'d';'e';'f';'g';'h';'i'};
s=1;
[n,m]=size(IM_c);
for i=1:n
    for j=1:m
        plot([IM_c(i,j),IM_c(i,j)],[0,1],'--b','linewidth',1.5);
        txt=text(IM_c(i,j)-0.06,0.1,['\itIM_',IM_subscript{s},' \rm= ',num2str(round(IM_c(i,j)*100)/100)]);
        set(txt,'Rotation',90);
        s=s+1;
    end
end
axis([0 IM_end 0 1.0]);
xlabel('Intensity Measure (IM)');
ylabel('Fragility Probability');
title('Original Fragility Curves and IM Level Candidates');
box on;

%% Output of the IM Level Candidates
disp('IM level candidates:')
y = IM_c

end

