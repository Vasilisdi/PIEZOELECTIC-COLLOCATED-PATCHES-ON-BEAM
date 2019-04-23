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
Length=0.01;         %length
epsilon33=1800*8.85*1e-12;       %clamped dielectric constant  
Cp=Length*Wp*epsilon33/tpax;

khatnum=E*I*3*Epatch*(((h/2)+tpax)^2-(h/2)^2)*d31;
khatden=4*tpax*(Epatch*((h/2+tpax)^3-(h/2)^3)+E*(h/2)^3);
khat=khatnum/khatden;

ga=khat/(p*Asurface);gs=(-d31*Epatch*Wp*((h/2)+tpax))/Cp;



zeta=1e-12*ones(8,1);
zeta2=[0.01;0.009;0.009;0.009;0.009;0.009;0.009;0.009];

%coef10e=1e+3;
coef10e=1e5;le2=1e+2;

n=30;                                      %seperation , #of elements

force=.1; %Nt
%lQrcon2=1e+20;lR=1e+5;         %ayto apodykneiei oti synadei me tis alles methodoys
lQrcon2=1e2;lR=1e+7;

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
temp=2;temp2=1;temp3=1;holdxs=zeros(numofpatches,2);bi=zeros(numofpatches,size(Mr,1));
for i=1:numofpatches
b(temp,temp2)=-1;b(temp+2,temp2)=1;
bi(temp2,4*i-2)=-1;bi(temp2,4*i)=1;
holdxs(i,[1 2])=[(ltot/n)*temp3 (ltot/n)*(temp3+1)];
patches(i)=(temp/2)+1;
temp=temp+4;temp2=temp2+1;
temp3=temp3+2;
end
%%
%show patch position
figure;lin=0;lfin=ltot/n;
for i=1:n
plot(linspace(lin,lfin,10),zeros(1,length(linspace(lin,lfin,10))),'b');hold on
if mod(i,2)==0 && i<2*numofpatches+1; plot(linspace(lin,lfin,10),zeros(1,length(linspace(lin,lfin,10))),'r','linewidth',2);hold on;end
lin=lfin;lfin=lfin+ltot/n; 
end;xlim([0 ltot]);set(gca,'YTick',[]);xlabel('x[m]')


%%

SQ=sqrt(valsqr);Z=zeros(size(SQ,1),size(SQ,2));
for i=1:2*n
if i==1;Z(i,i)=SQ(i,i)*2*zeta(1,1);end
if i==2;Z(i,i)=SQ(i,i)*2*zeta(2,1);end
if i~=1 && i~=2;Z(i,i)=SQ(i,i)*2*zeta(3,1);end
end
Cr=Mr*vec*Z*vec'*Mr;


A=[zeros(2*n,2*n),eye(2*n);-inv(Mr)*Kr,-inv(Mr)*Cr];
B=[zeros(2*n,2*n);inv(Mr)]*b;
C=[gs*bi,zeros(size(bi))];Dss=0;
Bw=[zeros(2*n,2*n);inv(Mr)]*bw;

sysi=ss(A,[B Bw],C,0);sysi.OutputName='usens.';



%%


white_noise(1,:)=force*zeros(1,1e5); 
white_noise(1,1:40000)=force*ones(1,40000);
tnoise=linspace(0,5,size(white_noise,2));


figure;lsim(sysi(:,9),white_noise,tnoise);title('no actuators included');

%%

load('Preumondmodel.mat')
Q=1e-2;R=1e7;
clear lqr;lqr=lqr(Model.a,Model.b(:,1:8),Q*eye(size(Model.a,1)),R*eye(size(Model.b(:,1:8),2)));

%define Lp observer matrix
% temp=1;
% poles2=-1*1e5*(1:size(Model.a,1)/2);poles=zeros(1,2*size(poles2,2));
% for i=1:length(poles2)
%    poles(1,[temp,temp+1])=[poles2(1,i)+4e4i;poles2(1,i)-4e4i];
%    temp=temp+2;
% end
poles=-1e6*(1:size(Model.a,1));
Lp=place(Model.a',Model.c',poles)';

Bm=Model.b(:,1:end-1);Bwm=Model.b(:,end);

noi=0.3;
%%
a=[Model.a-Lp*Model.c-Bm*lqr Lp*C;-B*lqr A];                  %1:16 Model , 17:32 Plant
b=[Lp;zeros(120,8)]; %system,measurement,control actions,force
c=[Model.c,zeros(size(C));zeros(size(Model.c)),C];              
d=0;
sys=ss(a,b,c,d);
figure;pzmap(sys)