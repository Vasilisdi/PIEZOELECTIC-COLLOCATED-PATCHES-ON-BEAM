clear all;clc;close all;   %#ok<CLSCR>
%%

%important  parameters

%my beam
%p=7880;         %steel
%p=7880;         %steel
%p=1800 ;         %epoxy resin/fiber glass beam
p=2770 ;         %aluminium
Width=0.02;h=0.003;
Asurface=Width*h;
I=(1/12)*Width*(h^3);
%E=2.1*(10^11);       %steel
%E=60*(10^9);          %epoxy regin/fiber glass beam
E=70*(10^9);          %alum
ltot=0.3;
m=p*Asurface*ltot;


%my patch
Epatch=50*1e9;             %piezoelectric patch Young mod piezo 50e9.
d31=-150*1e-12;               %piezoelectric patch
Wp=0.015;             %width
tpax=0.00065;          %thickness
e31=-7.5;
Length=0.012;         %length
epsilon33=1800*8.85*1e-12;       %clamped dielectric constant  
Cp=Length*Wp*epsilon33/tpax;

khatnum=E*I*3*Epatch*(((h/2)+tpax)^2-(h/2)^2)*d31;
khatden=4*tpax*(Epatch*((h/2+tpax)^3-(h/2)^3)+E*(h/2)^3);
khat=khatnum/khatden;

ga=khat/(p*Asurface);gs=(-d31*Epatch*Wp*((h/2)+tpax))/Cp;



zeta=5e-10*ones(8,1);
zeta2=0.1*[0.005;0.005;0.005;0.005;0.005;0.005;0.005;0.005];

%coef10e=1e+3;
coef10e=1e5;le2=1e+2;

n=25;                                      %seperation , #of elements

force=.1; %Nt
%lQrcon2=1e+20;lR=1e+5;         %ayto apodykneiei oti synadei me tis alles methodoys
lQrcon2=1e2;lR=1e7;         %ayto apodykneiei oti synadei me tis alles methodoys
%lQrcon2=1e5;lR=1e3;         %ayto apodykneiei oti synadei me tis alles methodoys

%%
%FEMs
%consists of #n FEMs


%determine the number of the elements seperation

l=ltot/n;
ndofs=2*(n+1);                             %  # of dofs
Mtot=zeros(ndofs,ndofs);Ktot=zeros(ndofs,ndofs);
temp=0;kE=(E*I)/(l^3);mE=(1/420)*(p*Asurface*l);

for i=1:n                                 %Mtot & Ktot
    K=kE*([12 6*l -12 6*l ; 6*l 4*(l^2) -6*l 2*(l^2) ; -12 -6*l 12 -6*l ; 6*l 2*(l^2) -6*l 4*(l^2)]);
    M=mE*([156 22*l 54 -13*l ; 22*l 4*(l^2) 13*l -3*(l^2) ; 54 13*l 156 -22*l ; -13*l -3*(l^2) -22*l 4*(l^2)]);
    Mtot(temp+1:temp+4,temp+1:temp+4)=Mtot(temp+1:temp+4,temp+1:temp+4)+M;
    Ktot(temp+1:temp+4,temp+1:temp+4)=Ktot(temp+1:temp+4,temp+1:temp+4)+K;
    temp=temp+2;
end

%%
%contraints determination for fems
%eigenvalues & eigenvectors

Kr=Ktot(3:size(Ktot,1),3:size(Ktot,2));Mr=Mtot(3:size(Ktot,1),3:size(Ktot,2));
[vec,valsqr]=eig(Kr,Mr);Eigenvalues=sqrt(diag(valsqr));

%%
%modes
le=linspace(0,ltot/n,15);               %x-axis position
Wh=zeros(n,length(le),size(vec,2));    %modes matrix
vec1=[zeros(2,size(vec,1));vec];
for z=1:(size(vec,2))                  %for every mode 
      temp=0;
 for j=1:n                             %for every part

  for i=1:length(le)                   %for part length(plot)
  
   phi1=1-3*(le(i)/(ltot/n))^2+2*(le(i)/(ltot/n))^3;
   phi2=(ltot/n)*((le(i)/(ltot/n))-2*(le(i)/(ltot/n))^2+(le(i)/(ltot/n))^3);
   phi3=3*(le(i)/(ltot/n))^2-2*(le(i)/(ltot/n))^3;
   phi4=(ltot/n)*(-(le(i)/(ltot/n))^2+(le(i)/(ltot/n))^3);
   Wh(j,i,z)=phi1*vec1(1+temp,z)+phi2*vec1(2+temp,z)+phi3*vec1(3+temp,z)+phi4*vec1(4+temp,z);
  
  end 
      temp=temp+2;
 end
end

%%
%modeshapes WH 

WH=zeros(size(Wh,3),size(Wh,2)*size(Wh,1));

for i=1:size(Wh,3)
    temp=0;
   for j=1:size(Wh,1)
      WH(i,1+temp:size(Wh,2)+temp)=Wh(j,:,i);
      temp=temp+size(Wh,2) ;
   end
end
lee=[];constant=0;
for i=1:size(Wh,1)
    lee=[lee le+constant];                                %#ok<AGROW>
    constant=constant+(ltot/n);
end

figure;
for i=1:size(WH,1)
    subplot(size(WH,1)/2,2,i);plot(lee,WH(i,:));
    ylabel(['X_' num2str(i), '[m]']);xlabel(['x[m]  ' , '-  mode for eigenvalue  :  ',num2str(Eigenvalues(i))])
end
clf;close figure 1
%%
%Patches

temp=2;numofpatches=0;
for i=1:floor(n/2)
temp=temp+4;numofpatches=numofpatches+1;
end
if n>8; numofpatches=8; end
b=zeros(2*n,numofpatches);             %determine bs the refering to control actions on evr dof - patches implementation
bw=b(:,1);bw(size(bw,1)-1,1)=1;        %determine the vibration - force at the right edge

patches=zeros(1,numofpatches);
temp=2;temp2=1;temp3=1;holdxs=zeros(numofpatches,2);
for i=1:numofpatches
b(temp,temp2)=-1;b(temp+2,temp2)=1;
holdxs(i,[1 2])=[(ltot/n)*temp3 (ltot/n)*(temp3+1)];
patches(i)=(temp/2)+1;
temp=temp+4;temp2=temp2+1;
temp3=temp3+2;
end
bi=b;
b=b*ga;
SQ=sqrt(valsqr);Z=zeros(size(SQ,1),size(SQ,2));
for i=1:2*n
if i==1;Z(i,i)=SQ(i,i)*2*zeta(1,1);end
if i==2;Z(i,i)=SQ(i,i)*2*zeta(2,1);end
if i~=1 && i~=2;Z(i,i)=SQ(i,i)*2*zeta(3,1);end
end
Cr=Mr*vec*Z*vec'*Mr;

%%
%show patch position
figure;lin=0;lfin=ltot/n;
for i=1:n
plot(linspace(lin,lfin,10),zeros(1,length(linspace(lin,lfin,10))),'b');hold on
if mod(i,2)==0 && i<2*numofpatches+1; plot(linspace(lin,lfin,10),zeros(1,length(linspace(lin,lfin,10))),'r','linewidth',2);hold on;end
lin=lfin;lfin=lfin+ltot/n; 
end;xlim([0 ltot]);set(gca,'YTick',[]);xlabel('x[m]')



%%


Mhold=eye(length(valsqr));                 %Mhold=vec'*Mr*vec check
SQ=sqrt(valsqr);Z=zeros(size(SQ,1),size(SQ,2));
for i=1:length(valsqr)
if i==1;Z(i,i)=SQ(i,i)*2*zeta(1,1);end
if i~=1;Z(i,i)=SQ(i,i)*2*zeta(2,1);end
end
K=zeros(length(valsqr),length(valsqr));
for i=1:length(valsqr)
    %K(i,i)=em(i)*(omega(i)^2);                  %preumond - natsiavas
    K(i,i)=valsqr(i,i);
end



%%
%input Vin and output u0
order=3;
Polyd=zeros(numofpatches,order,size(vec,2));holdit=zeros(size(vec,2),numofpatches);holditu0=zeros(size(vec,2),numofpatches);
for z=1:(size(vec,2))                  %for every mode 
 for j=1:numofpatches                  %for every part
   
      Polyd(j,:,z)=polyder(polyfit(le,Wh(2*j,:,z),order));

      holdit(z,j)=ga*(polyval(Polyd(j,:,z),le(end))-polyval(Polyd(j,:,z),le(1)));
      holditu0(z,j)=(polyval(Polyd(j,:,z),le(end))-polyval(Polyd(j,:,z),le(1)));  %8x60
      
 end
end

holdit2u0=(gs)*holditu0';

%check if its same O*vec
Mnew=Mhold;
Znew=Z;
Knew=K;

Anew3=[zeros(numofpatches,numofpatches)  , eye(numofpatches) ;  ...
    -(Mnew(1:8,1:8)/holdit2u0(1:8,1:8))\(Knew(1:8,1:8)/holdit2u0(1:8,1:8)) , ...
    -(Mnew(1:8,1:8)/holdit2u0(1:8,1:8))\(Znew(1:8,1:8)/holdit2u0(1:8,1:8))];
vecb=vec'*b;
Bnew3=[zeros(numofpatches,numofpatches) ; (Mnew(1:8,1:8)/holdit2u0(1:8,1:8))\vecb(1:8,1:8)];
vecbw=vec'*bw;
Bwnew3=[zeros(numofpatches,1) ; (Mnew(1:8,1:8)/holdit2u0(1:8,1:8))\vecbw(1:8,1)];
Cnew3=eye(2*numofpatches);Dnew3=0;

sys2=ss(Anew3,[Bnew3 Bwnew3],Cnew3,Dnew3);sys3new.OutputName='u out';
sys2.OutputName='u out';


Qlqr2=eye(2*numofpatches);Qlqr2=lQrcon2*Qlqr2;
%Qlqr2=zeros(2*numofpatches,2*numofpatches);for i=1:numofpatches;Qlqr2(i,i)=1;end;Qlqr2=lQrcon2*Qlqr2;
R2=diag(lR*ones(1,numofpatches));N2=zeros(size(Bnew3,1),length(R2));

[Kcontroller2,~,e]=lqr(Anew3,Bnew3,Qlqr2,R2,N2);


%check cotrollability
% check=[Qlqr2 N2;N2' R2];positivity=eig(check);
% inv(check);
% chol(check);
% 


% clear Pcontrolability 
% Pcontrollability=ctrb(Anew3,Bnew3);
% unco2=length(Anew3)-rank(Pcontrollability);

%any piezopatch applied 

%piezoelectric tansducers included
sys2controlled=ss(Anew3-Bnew3*Kcontroller2,Bwnew3,Cnew3,0);
sys2controlled.OutputName='u out';

white_noise(1,:)=force*zeros(1,1e5); 
white_noise(1,1:10000)=force*ones(1,10000);
tnoise=linspace(0,5,size(white_noise,2));


figure;subplot(1,2,1);lsim(sys2(1:numofpatches,numofpatches+1),white_noise,tnoise);title('no actuators included');
   subplot(1,2,2);lsim(sys2controlled(1:numofpatches,1),white_noise,tnoise);title('transducers included - LQR')

figure;subplot(1,2,1);pzmap(sys2*force);title('initial system');grid on;grid on
subplot(1,2,2);pzmap(sys2controlled*force);title('c.l.sys2. - lqr');grid on
%figure;pzmap(sys2*force);title('initial system');grid on;xlim([-1 1e-8])
          
figure;lsim(-Kcontroller2*sys2controlled,white_noise,tnoise);title('state feedback for lqr ')





%%

%MPC for both models




%MPC for both models


%%
%MPC application for both models
coef=0.3;
SampleTime=0.025;

%first application
%application to initial system with time constants of modes outputs 

%application to initial system with voltage of sensors outputs 
sys3=ss(sys2.a,sys2.b,eye(size(sys2.a,1)),0);
clear Am Bm Cm constantofR R Q Np Nc u t y
sysd2=c2d(sys3(1:8,:),SampleTime,'zoh');
[Ad2,Bd2,Cd2,Dd2,Kd2]=ssdata(sysd2);
Np=80;  %constant number
Nc=23;  %constant number
constantofR=.04;constantofQ=1;
%constantofR=0;constantofQ=1;
Am=Ad2;Bm=Bd2(:,1:numofpatches);Bmw=Bd2(:,numofpatches+1);Cm=Cd2;             
R=eye(Nc*size(Bm,2))*constantofR;Q=constantofQ*eye(Np*size(Cm,1));  %Q(NpxNp) και R(NcxNc) , cost of deltaXs and deltaUs respectively.
[ sysMPCd2 , Kmpc2] = SyscratorWITHMPC( Np , Nc , R , Q , Am , Bm , Cm , SampleTime);


sysMPCd2=ss(sysMPCd2.a,[zeros(16,8),Bmw;eye(8),Cm*Bmw],sysMPCd2.c,0,SampleTime);
%find Ysp, input signal
%tfin=10;
%[Ysp,t]=step(sys3new(:,numofpatches+1)*force,0:SampleTime:tfin);




tfin=100;

tn=0:SampleTime:tfin;len=length(tn);forcempc=zeros(size(tn));
force2mpc=force*(sin(2*tn));
for i=2:length(tn)
forcempc(i)=force*(sin(2*tn(i))-sin(2*tn(i-1)));
end;forcempc(1)=0;
whitensys=coef*2*5e-3*randn(8,len);
whitenoisemeas=coef*1e-3*randn(8,len);



figure;[yr,tr]=lsim(sysMPCd2(17:24,:),[whitensys;forcempc],tn);
for i=1:size(yr,2);   subplot(size(yr,2),1,i); stairs(tr,yr(:,i));if i~=8; set(gca,'XTick',[]); end ;if i==1;title('state vector');end;end

figure
for i=1:size(yr,2);  subplot(size(yr,2),1,i); stairs(tr,yr(:,i)+whitenoisemeas(i,:)');
    if i~=8; set(gca,'XTick',[]); end ; 
    if i==1;title('measurements with noise');end;
end

figure;pzmap(sysMPCd2);grid on; %xlim([-1 1])

%feedbackcontrol

figure;[uk,tk]=lsim(-Kmpc2*sysMPCd2,[whitensys;forcempc],tn);
for i=1:size(uk,2); subplot(size(uk,2),1,i); stairs(tk,uk(:,i));
    if i~=8; set(gca,'XTick',[]); end ; 
    if i==1;title('input vector');end;
    ylabel(['V_{in} _' ,num2str(temp)]); 
grid on
end
%%
%observer
% load('Modelmpc.mat')

A=[Am zeros(size(Am,1),size(Cm,1));Cm*Am eye(size(Cm,1))];
B=[Bm ; Cm*Bm];
E=sysMPCd2.b;
 C=eye(size(A,1));
 %Lp=lqr(A',C',eye(size(A,1)),1000*eye(size(C,1)))';Lp=1e-3*Lp(:,17:17+7);
%Lp=place(A',C',0.254*pole(sysMPCd2))';Lp=Lp(:,17:17+7);

C=[zeros(8,16),eye(8)];
Lp=place(A',C',0.95*pole(sysMPCd2))';
%feedback controlled with estimator
a=[A -B*Kmpc2;Lp*C A-Lp*C-B*Kmpc2];                  %1:24 Model , 25:48 Plant
b=[E,zeros(size(Lp));zeros(size(E)),Lp];             %system noise, force and measurement-output noise respectively
d=0;
sys=ss(a,b,eye(size(a,1)),d,SampleTime);
dist=[whitensys;forcempc;whitenoisemeas];



[y2,t2]=lsim(sys(41:48,:),dist',tn);                   %compensation with estimator

figure;temp=1;
for i=1:8
subplot(8,1,temp);stairs(t2,y2(:,i));
if temp~=8; set(gca,'XTick',[]); end ; 
ylabel(['u_{out} _' ,num2str(temp)]);
if temp==1; title('state vector'); end; temp=temp+1;
grid on
end;

figure;temp=1;
for i=1:8
subplot(8,1,temp);stairs(t2,y2(:,i)+whitenoisemeas(i,:)');
if temp~=8; set(gca,'XTick',[]); end ; 
ylabel(['u_{out} _' ,num2str(temp)]);
if temp==1; title('measuments including noise'); end; temp=temp+1;
grid on
end;

figure;pzmap(sys);grid on;  %xlim([-1 1])

[u2,t2]=lsim(-Kmpc2*sys(1:24,:),dist',tn); 

figure;
for i=1:size(u2,2)
subplot(size(u2,2),1,i);stairs(t2,u2(:,i)')
if i~=8; set(gca,'XTick',[]); end ; 
ylabel(['V_{in} _' ,num2str(temp)]);
if i==1; title('input vector'); end; 
grid on
end;
%%
figure
temp=1;
for i=1:8
    subplot(8,2,temp);stairs(tr,yr(:,i))
    if temp~=15; set(gca,'XTick',[]); end ;
    if temp==1;title('state vector only mpc');end
    subplot(8,2,temp+1);stairs(t2,y2(:,i))
    if temp~=15; set(gca,'XTick',[]); end ;
    if temp==1;title('state vector including observer');end
   temp=temp+2; 
end

%%
%noise

figure;temp=1;
for i=1:8;
    subplot(8,2,temp);plot(tn,whitensys(i,:));if i~=8; set(gca,'XTick',[]); end ;grid on
    if i==1; title('system white noise'); end;if i==8; xlabel('t[s]'); end ;
    ylabel(['u_{out} _',num2str(i) ,'[V]'])
    subplot(8,2,temp+1);plot(tn,whitenoisemeas(i,:));if i~=8; set(gca,'XTick',[]); end ;grid on
    if i==1; title('measurement white noise'); end;if i==8; xlabel('t[s]'); end ;
    ylabel(['u_{out} _',num2str(i) ,'[V]'])

    temp=temp+2;
end
%%

%sysreal view


A=[zeros(numofpatches,numofpatches)  , eye(numofpatches) ; -Mnew(1:8,1:8)\Knew(1:8,1:8) , -Mnew(1:8,1:8)\Znew(1:8,1:8)];
B=[zeros(numofpatches,numofpatches) ; Mnew(1:8,1:8)\vecb(1:8,1:8)];    %for state space control
Bdw=[zeros(numofpatches,1) ; Mnew(1:8,1:8)\vecbw(1:8,1)]; 
C=[eye(numofpatches),zeros(numofpatches,numofpatches)];
D=0;

sysr=ss(A,[B,Bdw],C,D);

[uact,tact]=lsim(-Kcontroller2*sys2controlled,white_noise,tnoise);

[xs,ts]=lsim(sysr,[uact';white_noise],tnoise);

figure;plot(ts,xs*WH(1:numofpatches,end));title('right end of beam')


%%

clear Am Bm Cm 
sysd2=c2d(sysr(1:8,:),SampleTime,'zoh');
[Ad2,Bd2,Cd2,Dd2,Kd2]=ssdata(sysd2);

Am=Ad2;Bm=Bd2(:,1:numofpatches);Bmw=Bd2(:,numofpatches+1);Cm=Cd2; 
A=[Am zeros(size(Am,1),size(Cm,1));Cm*Am eye(size(Cm,1))];
B=[Bm ; Cm*Bm];
E=[zeros(16,8),Bmw;Cm*Bm,Cm*Bmw];
sysr2=ss(A,[E,B],[zeros(8,16),eye(8)],0,SampleTime);

[xs,ts]=lsim(sysr2,[whitensys;forcempc;uk'],tn);

figure;stairs(ts,xs*WH(1:8,end));hold on;title('right end of beam')
xlabel('t[s]');ylabel('x[m]');xlim([0 10])

[xs,ts]=lsim(sysr2,[whitensys;forcempc;u2'],tn);
stairs(ts,0.5*xs*WH(1:8,end))
legend('mpc','compensator')

figure;

[xs,ts]=lsim(sysr2,[whitensys;forcempc;0*u2'],tn);
stairs(ts,0.5*xs*WH(1:8,end))
xlabel('t[s]');ylabel('x[m]');xlim([0 10])
