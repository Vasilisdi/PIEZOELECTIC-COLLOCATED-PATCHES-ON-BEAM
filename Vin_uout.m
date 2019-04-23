close all;clear all;clc;         %#ok<CLSCR>
%%
load('lqr.mat')            %lqr from FEM model
load('Preumondplant.mat')  %Plant
load('Preumondmodel.mat')  %Model
load('Preumondplant2.mat')
load('lqr2.mat')
load('lqr3.mat')           %for plant



lqr=lqr3;



%Lp criterion is below 

% temp=1;
% poles2=-5*1e5*linspace(2,6,size(Model.a,1)/2);poles=zeros(1,2*size(poles2,2));
% %poles2=-1*1e7*(1:size(Model.a,1)/2);poles=zeros(1,2*size(poles2,2));
% %poles2=-1*1e6*(1:size(Model.a,1)/2);poles=zeros(1,2*size(poles2,2));
% for i=1:length(poles2)
%    poles(1,[temp,temp+1])=[poles2(1,i)+0*1e4i;poles2(1,i)-0*1e4i];
%    %poles(1,[temp,temp+1])=[poles2(1,i)+0*1e3i;poles2(1,i)-0*1e3i];
%    temp=temp+2;
% end


poles=-.2*1e7*linspace(2,6,size(Model.a,1));
Lp=place(Model.a',Model.c',(poles))';

Bp=Plant.b(:,1:end-1);Bwp=Plant.b(:,end);
Bm=Model.b(:,1:end-1);Bwm=Model.b(:,end);

noi=0.06;
force1=.1;       %[Nt]   0.1Nt

tsim=5;
%%
a=[Model.a-Lp*Model.c-Bm*lqr Lp*Plant.c;-Bp*lqr Plant.a];                  %1:16 Model , 17:32 Plant
b=[zeros(size(Plant.c')),Lp,Bm,zeros(size(Bwm));Plant.c',zeros(size(Lp)),Bp,Bwp]; %system,measurement,control actions,force
c=[Model.c,zeros(size(Model.c));zeros(size(Plant.c)),Plant.c];              
d=0;

sys=ss(a,b,c,d);
%sys2=ss(a,b,eye(size(a,1)),d);
sys2=ss(a,b,eye(size(a,1)),d);

%%
%disturbances d white noise added to output
%len=3*1e4;tsim=1;             %define dist.
len=3*1e4;             %define dist.
white_noise1=zeros(size(Model.c,1),len);white_noise2=zeros(size(Model.c,1),len);white_noise3=zeros(size(Model.c,1),len);
for i=1:size(Model.c,1)
   white_noise1(i,:)=noi*0.005*randn(1,size(white_noise1,2));       %on ds system noise
   white_noise2(i,:)=noi*0.01*0*randn(1,size(white_noise2,2));         %on dy measurement noise
   white_noise3(i,:)=noi*0.01*randn(1,size(white_noise3,2));         %on du cotrol action noise
end
tnoise=linspace(0,tsim,size(white_noise1,2));
z=zeros(size(white_noise1));

%disturbance w if needed, interpretation : force at edge
%force added to input w 

w0=zeros(size(tnoise));
for i=1:floor(size(w0,2)/10)
w0(1,i)=force1;
end


%%
%the part below proves that the estimator application is not convenient for
%such a problem. in case force is involved


[y20,t20]=lsim(sys,[white_noise1;white_noise2;white_noise3;w0],tnoise);                   %compensation with estimator
%for control actions
[yu20,tu20]=lsim(sys2,[white_noise1;white_noise2;white_noise3;w0],tnoise);                %compensation with estimator
 


%the control actions
%control actions, estimator included
u20=zeros(size(Bm,2),length(tu20));
for i=1:length(tu20)
u20(:,i)=-lqr*(yu20(i,1:size(Model.a))')+white_noise3(:,i);   
end 


%%
%simulation noise and force



figure;temp=1;
for i=1:size(Model.c,1);
    subplot(size(Model.c,1),3,temp);stairs(tnoise,white_noise1(i,:));if i~=size(Model.c,1); set(gca,'XTick',[]); end ;grid on
    if i==1; title('system white noise'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
    subplot(size(Model.c,1),3,temp+1);stairs(tnoise,white_noise2(i,:));if i~=size(Model.c,1); set(gca,'XTick',[]); end ;grid on
    if i==1; title('measurement white noise'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
    subplot(size(Model.c,1),3,temp+2);stairs(tnoise,white_noise3(i,:));if i~=size(Model.c,1); set(gca,'XTick',[]); end ;grid on
    if i==1; title('control action white noise'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;
    ylabel(['u_{out} _',num2str(i) ,'[V]'])

    temp=temp+3;
end




figure;temp=1;
for i=1:size(Model.c,1);
    subplot(size(Model.c,1),1,temp);plot(t20,y20(:,i+size(Model.c,1)));if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('state vars.'); end; if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
     
    temp=temp+1;
end 

figure;temp=1;
for i=1:size(Model.c,1);
    subplot(size(Model.c,1),1,temp);plot(t20,y20(:,i+size(Model.c,1))+white_noise2(i,:)');if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('measurements (including noise)'); end; if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
     
    temp=temp+1;
end 


%simulate control actions
figure;temp=1;
for i=1:size(u20,1)
    
    
    subplot(size(u20,1),1,temp);plot(tu20,u20(i,:));                         %u2 control actions with estimator
    if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('input vector'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['V_{in} _',num2str(i) ,'[V]'])
    
    temp=temp+1;
end


%%
%poles of system

figure;
pzmap(sys);title('compensationed system');grid on
figure;pzmap(ss(Plant.a-Bp*lqr,Bp,eye(size(Plant.a)),0));grid on;title('only lqr');

%%

a=Plant.a;                  %1:16 Model , 17:32 Plant
%b=[Lp,Bwm;zeros(size(Lp)),Bwp];
b=[Plant.c',Bp,Bwp];        %system,control actions,force
c=Plant.c;
d=[eye(size(Plant.c,1)),zeros(size(Plant.c,1),size([Bm,Bwm],2))];
%d=[zeros(size(Plant.c,1),size(Plant.c,1)),zeros(size(Plant.c,1),1);eye(size(Plant.c,1)),zeros(size(Plant.c,1),1)];
sysMod=ss(a-Bp*lqr,b,c,d);

figure;[y,t]=lsim(sysMod,[white_noise1;white_noise3;w0],tnoise);
temp=1;
for i=1:size(Model.c,1);
    subplot(size(Model.c,1),1,temp);stairs(t20,y20(:,i+size(Model.c,1)),'b','LineWidth',0.5);hold on;stairs(t,y(:,i),'r','LineWidth',0.1)
    if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('comparison'); end; if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
    temp=temp+1;
    if i==1;legend('compensationed plant','lqr');end
    xlim([4.9 5])
end

%the control actions
%control actions, estimator included
sysMod=ss(sysMod.a,sysMod.b,eye(size(sysMod.a,1)),0);
[u201,t201]=lsim(-lqr*sysMod,[white_noise1;white_noise3;w0],tnoise);



figure;temp=1;
for i=1:size(u201,2)
    
    
    subplot(size(u201,2),1,temp);stairs(t201,u201(:,i));                         %u2 control actions with estimator
    if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('input vector'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['V_{in} _',num2str(i) ,'[V]'])
    
    temp=temp+1;
end

%%


load('Wnew.mat')

sysr=Plant2;

sysr=ss(sysr.a,[Plant.c',sysr.b],sysr.c,0); %system noise,control actions,force
[xs,ts]=lsim(sysr,[white_noise1',u20',w0']',tnoise);

figure;plot(ts,0.9*xs*Wnew(:,end),'LineWidth',1.2);title('right end of beam');hold on;xlim([0 2])

[xs,ts]=lsim(sysr,[white_noise1',u201,w0']',tnoise);
plot(ts,xs*Wnew(:,end));legend('compensator','only lqr')
xlabel('t[s]');ylabel('x[m]');

figure;[xs,ts]=lsim(sysr,[white_noise1',0*u201,w0']',tnoise);
plot(ts,xs*Wnew(:,end));title('no active vibration control');xlabel('t[s]');ylabel('x[m]');

stoop
%%


%DISCRETE SYSTEM




%%
T=0.0025;


sys=c2d(sys,T);
%sys2=ss(a,b,eye(size(a,1)),d);
sys2=c2d(sys2,T);

%%
%disturbances d white noise added to output
%len=3*1e4;tsim=1;             %define dist.
t=0:T:tsim;len=length(t);
white_noise1=zeros(size(Model.c,1),len);white_noise2=zeros(size(Model.c,1),len);white_noise3=zeros(size(Model.c,1),len);
for i=1:size(Model.c,1)
   white_noise1(i,:)=noi*0.005*randn(1,size(white_noise1,2));       %on ds system noise
   white_noise2(i,:)=noi*0*randn(1,size(white_noise2,2));         %on dy measurement noise
   white_noise3(i,:)=noi*0.01*randn(1,size(white_noise3,2));         %on du cotrol action noise
end
tnoise=t;
z=zeros(size(white_noise1));

%disturbance w if needed, interpretation : force at edge
%force added to input w 

w0=zeros(size(tnoise));
for i=1:0.5/T
w0(1,i)=force1;
end


%%
%the part below proves that the estimator application is not convenient for
%such a problem. in case force is involved


[y20,t20]=lsim(sys,[white_noise1;white_noise2;white_noise3;w0],tnoise);                   %compensation with estimator
%for control actions
[yu20,tu20]=lsim(sys2,[white_noise1;white_noise2;white_noise3;w0],tnoise);                %compensation with estimator
 


%the control actions
%control actions, estimator included
u20=zeros(size(Bm,2),length(tu20));
for i=1:length(tu20)
u20(:,i)=-lqr*(yu20(i,1:size(Model.a))')+white_noise3(:,i);   
end 


%%
%simulation noise and force



figure;temp=1;
for i=1:size(Model.c,1);
    subplot(size(Model.c,1),3,temp);stairs(tnoise,white_noise1(i,:));if i~=size(Model.c,1); set(gca,'XTick',[]); end ;grid on
    if i==1; title('system white noise'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
    subplot(size(Model.c,1),3,temp+1);stairs(tnoise,white_noise2(i,:));if i~=size(Model.c,1); set(gca,'XTick',[]); end ;grid on
    if i==1; title('measurement white noise'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
    subplot(size(Model.c,1),3,temp+2);stairs(tnoise,white_noise3(i,:));if i~=size(Model.c,1); set(gca,'XTick',[]); end ;grid on
    if i==1; title('control action white noise'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;
    ylabel(['u_{out} _',num2str(i) ,'[V]'])

    temp=temp+3;
end




figure;temp=1;
for i=1:size(Model.c,1);
    subplot(size(Model.c,1),1,temp);stairs(t20,y20(:,i+size(Model.c,1)));if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('state vars.'); end; if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
     
    temp=temp+1;
end 

figure;temp=1;
for i=1:size(Model.c,1);
    subplot(size(Model.c,1),1,temp);stairs(t20,y20(:,i+size(Model.c,1))+white_noise2(i,:)');if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('measurements (including noise)'); end; if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
     
    temp=temp+1;
end 


%simulate control actions
figure;temp=1;
for i=1:size(u20,1)
    
    
    subplot(size(u20,1),1,temp);stairs(tu20,u20(i,:));                         %u2 control actions with estimator
    if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('input vector'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['V_{in} _',num2str(i) ,'[V]'])
    
    temp=temp+1;
end


%%
%poles of system

figure;
pzmap(sys);title('compensationed system');grid on
sysMod=c2d(sysMod,T);
figure;pzmap(sysMod);grid on;title('only lqr');
sysMod=ss(sysMod.a,sysMod.b,eye(size(sysMod.c,2)),0,T);
%%




figure;[y,t]=lsim(sysMod(1:8,:),[white_noise1;white_noise3;w0],tnoise);
temp=1;
for i=1:size(Model.c,1);
    subplot(size(Model.c,1),1,temp);stairs(t20,y20(:,i+size(Model.c,1)),'b','LineWidth',0.5);hold on;stairs(t,y(:,i),'r','LineWidth',0.1)
    if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('comparison'); end; if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
    temp=temp+1;
    if i==1;legend('compensationed plant','lqr');end
    xlim([0.4 .5])
end


%the control actions
%control actions, estimator included
u201=zeros(size(Bm,2),length(t));[y,t]=lsim(sysMod,[white_noise1;white_noise3;w0],tnoise);
for i=1:length(t)
u201(:,i)=-lqr*(y(i,1:size(Model.a))')+white_noise3(:,i);   
end 


figure;temp=1;
for i=1:size(u201,1)
    
    
    subplot(size(u201,1),1,temp);stairs(t,u201(i,:));                         %u2 control actions with estimator
    if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
    if i==1; title('input vector'); end;if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
    ylabel(['V_{in} _',num2str(i) ,'[V]'])
    
    temp=temp+1;
end


% %%
% %comparison
% 
% figure;[y,t]=lsim(sysMod(1:8,:),[white_noise1;white_noise3;w0],tnoise);
% temp=1;
% for i=1:size(Model.c,1);
%     subplot(size(Model.c,1),1,temp);stairs(t20,y20(:,i+size(Model.c,1)),'b','LineWidth',0.5);hold on;stairs(t,y(:,i),'r','LineWidth',0.1)
%     if i~=size(Model.c,1); set(gca,'XTick',[]); end ;
%     if i==1; title('comparison'); end; if i==size(Model.c,1); xlabel('t[s]'); end ;grid on
%     ylabel(['u_{out} _',num2str(i) ,'[V]'])
%     temp=temp+1;
%     if i==1;legend('compensationed plant','lqr');end
%     xlim([9.5 9.9])
% end
% 

%%

sysr=Plant2;
sysr=c2d(ss(sysr.a,[Plant.c',sysr.b],sysr.c,0),T); %system noise,control actions,force
[xs,ts]=lsim(sysr,[white_noise1',u20',w0']',tnoise);

figure;stairs(ts,0.9*xs*Wnew(:,end),'LineWidth',1.2);title('right end of beam');hold on;

[xs,ts]=lsim(sysr,[white_noise1',u201',w0']',tnoise);
stairs(ts,xs*Wnew(:,end));legend('compensator','only lqr');xlim([0 2]);grid on
xlabel('t[s]');ylabel('x[m]');



[xs,ts]=lsim(sysr,[white_noise1',0*u201',w0']',tnoise);
figure;stairs(ts,xs*Wnew(:,end));title('right end of beam - including white noise');grid on
xlabel('t[s]');ylabel('x[m]');

