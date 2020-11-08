clc; clear;

% This is an example of fragility conversion (from Cloud to IDA for DS3)
% Damage State DS3 (extensive)
ds=3; IM_end=3.0;
im=linspace(0.001,IM_end,500);

% Target method (IDA)
eta_i =[0.5691;1.0896;1.6374;3.6147];
beta_i=[0.2800;0.3283;0.3978;0.4193];
pfi(:,ds)= normcdf((log(im/eta_i(ds)))/beta_i(ds));    

% Original method (Cloud)
eta_c =[0.5105;1.2304;2.1229;5.0267];
beta_c=[0.3886;0.3886;0.3886;0.3886];
pfc(:,ds)= normcdf((log(im/eta_c(ds)))/beta_c(ds)); 

% IM_candidate
y = IM_candidate(eta_c,beta_c,IM_end);

% Visual inspection by using Figure of IM candidate to select IM values
IM=[y(1),y(2),y(4),y(6)];

% Fragility points calculated by nonlinear time-history analysis
Aux_p1=[ 9/80,IM(2)];
Aux_p2=[74/80,IM(4)];

% Converted fragiltiy curve
im=linspace(0.001,IM_end,500);
[eta_u(ds),beta_u(ds)] = frag_conversion(Aux_p1(1,:),Aux_p2(1,:));
pfu(:,ds)= normcdf((log(im/eta_u(ds)))/beta_u(ds));

% Draw figure
figure
h1=plot(im,pfc(:,ds),'k','linewidth',1);
hold on
h2=plot(im,pfi(:,ds),'r','linewidth',2);
h3=plot(im,pfu(:,ds),'--b','linewidth',2);
h4=scatter(Aux_p1(1,2),Aux_p1(1,1),100,'b','linewidth',2);
h5=scatter(Aux_p2(1,2),Aux_p2(1,1),100,'r','linewidth',2);

legend([h1,h2,h3,h4,h5],'Cloud (original)','IDA (target)','Converted','IM1','IM2','Location','Best');
axis([0 IM_end 0 1.0]);
xlabel('Intensity Measure (IM)');
ylabel('Fragility Probability');
title('Original, target and converted Fragility Curves');
