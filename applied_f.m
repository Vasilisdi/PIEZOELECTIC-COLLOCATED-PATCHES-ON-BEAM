clear all; clc; close all; %#ok<CLSCR> 
%%
%FEM and modes
%beam FEM
%my beam
%p=7880;         %steel
p=2770 ;         %epoxy resin/fiber glass beam [kg/m3]
Width=0.02;h=0.003;
Asurface=Width*h;
I=(1/12)*Width*(h^3);
%E=2.1*(10^11);  %steel
E=70*(10^9);     %epoxy regin/fiber glass beam [N/m2]
ltot=.3;
m=p*Asurface*ltot;
%load
force=10;     %Nt
Length=0.01;
numofpatches=5;

n=10;                                      %seperation , #of elements
%determine the number of the elements seperation
%%


beta2=[1.875;4.694;7.855;10.996;14.137;17.279;(7-(1/2))*pi;(8-(1/2))*pi]./ltot;
beta=zeros(numofpatches,1);
if length(beta2)<numofpatches
    beta(1:length(beta2),1)=beta2;
    for i=length(beta2)+1:numofpatches
        beta(i,1)=(i-(1/2))*pi/ltot;
    end
end
if length(beta2)>=numofpatches
    beta=beta2(1:numofpatches,1);
end

%%

%determine zetas
%zeta=5e-10*ones(numofpatches,1);
zeta=[0.005;0.005;0.005;0.005;0.005;0.005;0.005;0.005];

%%
gama(:,1)=-(sin(beta*ltot)+sinh(beta*ltot))./(cos(beta*ltot)+cosh(beta*ltot));
x=linspace(0,ltot,301);W=zeros(size(beta,1),size(x,2));
omega=(sqrt(E*I/(p*Asurface))).*(beta.^2);
for i=1:size(x,2)
W(:,i)=sin(beta*x(i))-sinh(beta*x(i))+gama(:,1).*(cos(beta*x(i))-cosh(beta*x(i)));
end
%not the modes we have to find the an coef. of these modes  - no
%orthogonality conditions check
figure;
for i=1:size(W,1)
    if size(beta,1)==8 ; subplot(size(W,1)/2,2,i); end
    if size(beta,1)==3 ; subplot(size(W,1),1,i); end
    plot(x,W(i,:));ylabel(['X_' num2str(i), '[m]']);xlabel(['x[m] ' ,'-  mode for eigenvalue  :  ',num2str(omega(i)) ])
end
%%
clf;close figure 1


%%
%second mode pres. - to check the following seperation
figure;
for i=1:size(W,1)
    if size(beta,1)==8 ; subplot(size(W,1)/2,2,i); end
    if size(beta,1)==3 ; subplot(size(W,1),1,i); end
    plot(x,W(i,:));hold on;
end
%%
clf


%%

%seperate every mode into some parts

[nummodes,size2devide]=size(W);
orderp=4; %polyfit order
parts=10; %# of spited parts for seperation
forw=(size(x,2)-1)/parts; %# of used points
P=zeros(parts,orderp+1,nummodes);
xsp=zeros(parts,forw+1);temp=0;Wsp=zeros(parts,forw+1,nummodes);
for i=1:parts
xsp(i,:)=x(1,1+temp:forw+1+temp);
for j=1:nummodes
Wsp(i,:,j)=W(j,1+temp:forw+1+temp);
end
temp=temp+forw;
end

%expressed modes as polynomials
for i=1:nummodes
    for j=1:parts
    P(j,:,i)=polyfit(xsp(j,:),Wsp(j,:,i),orderp); 
    end
end
%check the seperation of modes with the above func. (in same plot)
for i=1:size(P,3)
    
    if size(beta,1)==8 ; subplot(size(P,3)/2,2,i);hold on; end
    if size(beta,1)==3 ; subplot(size(P,3),1,i);hold on; end
    
    for j=1:size(P,1)
    plot(xsp(j,:),polyval(P(j,:,i),xsp(j,:)));
    end
end
%%
clf;close figure 1


%%
%integrals and em. mass of modes
%check for the orthogonality conditions
em=zeros(nummodes,1);
for i=1:size(P,3)
    for j=1:size(P,1)
      em(i,1)=em(i,1)+diff(polyval(polyint(conv(P(j,:,i),P(j,:,i)),0),xsp(j,[1 forw+1])));
    end
end
em=em*p*Asurface;
modeeqs=sqrt(em);

%modes knowing an coefs. for every mode - an=1/sqrt(em)
Wnew=zeros(size(beta,1),size(x,2));
for i=1:size(x,2)
Wnew(:,i)=(1./modeeqs).*W(:,i);
end
%the modes knowing the an coefs  -  orthgonality checked
figure;
for i=1:size(Wnew,1)
    
    if size(beta,1)==8 ; subplot(size(Wnew,1)/2,2,i); end
    if size(beta,1)==3 ; subplot(size(Wnew,1),1,i); end
    
    plot(x,Wnew(i,:));ylabel(['X_' num2str(i), '[m]']);xlabel(['x[m] ' ,'-  mode for eigenvalue  :  ',num2str(omega(i)) ])
end
clf;close figure 1
%following the polynomial approach of the modes
figure;
for i=1:size(Wnew,1)
    if size(beta,1)==8 ; subplot(size(Wnew,1)/2,2,i); end
    if size(beta,1)==3 ; subplot(size(Wnew,1),1,i); end
    plot(x,Wnew(i,:));
end
%%
clf;
%%
%seperate in 8 parts in order to check the orthogonal conditions

xsp=zeros(parts,forw+1);temp=0;Wsp=zeros(parts,forw+1,nummodes);
for i=1:parts
xsp(i,:)=x(1,1+temp:forw+1+temp);
for j=1:nummodes
Wsp(i,:,j)=Wnew(j,1+temp:forw+1+temp);
end
temp=temp+forw;
end


%seperate modes in 8 parts
P=zeros(parts,orderp+1,nummodes);P2=zeros(parts,orderp-1,nummodes);
for i=1:nummodes
    for j=1:parts
    P(j,:,i)=polyfit(xsp(j,:),Wsp(j,:,i),orderp);  %polyfit for parts of modes
    P2(j,:,i)=polyder(polyder(P(j,:,i))); %2nd derivative of polyfit for parts of modes
    end
end
%check the seperation of modes with the above func.
for i=1:size(P,3)  
    if size(beta,1)==8 ; subplot(size(P,3)/2,2,i);hold on; end
    if size(beta,1)==3 ; subplot(size(P,3),1,i);hold on; end
    
    for j=1:size(P,1)
    plot(xsp(j,:),polyval(P(j,:,i),xsp(j,:)));
    end
end
%%
clf;close figure 1
%%

emnew=zeros(nummodes,1);
for i=1:size(P,3)
    for j=1:size(P,1)
      emnew(i,1)=emnew(i,1)+diff(polyval(polyint(conv(P(j,:,i),P(j,:,i)),0),xsp(j,[1 11])));
    end
end
emnew=emnew*p*Asurface;
%%
%equations of motion
M=zeros(nummodes,nummodes);Z=zeros(nummodes,nummodes);K=zeros(nummodes,nummodes);
for i=1:size(em,1)
    M(i,i)=1;
end
for i=1:size(em,1) 
    Z(i,i)=2*zeta(i)*omega(i);
end
for i=1:size(em,1)
    K(i,i)=(omega(i)^2);
end

%force

orderp=3;                             %polyfit order
F=zeros(nummodes,numofpatches);
PP=zeros(nummodes,orderp+1,numofpatches);
PForce=zeros(nummodes,1);
PForcen=zeros(nummodes,1);
temp=0;posdist=find(x(1,:)==Length-temp)-find(x(1,:)==0);
while isempty(posdist)==1
posdist=find(x(1,:)==Length-temp)-find(x(1,:)==0);
temp=0.001+temp;
end

figure;
for i=1:size(Wnew,1)
    
    if size(beta,1)==8 ; subplot(size(Wnew,1)/2,2,i); end
    if size(beta,1)==3 ; subplot(size(Wnew,1),1,i); end
    
    plot(x,Wnew(i,:));ylabel(['X_' num2str(i), '[m]']);xlabel(['x[m] ' ,'-  mode for eigenvalue  :  ',num2str(omega(i)) ])
end


figure;
for i=1:size(Wnew,1)
    subplot(size(Wnew,1),1,i);
    plot(x,Wnew(i,:));hold on
end

holdpatchpos=zeros(2,numofpatches);
for i=1:nummodes
    post=0;

    for j=1:numofpatches
        
    if j==1; post=posdist;end;if j~=1;post=2*posdist+post;end       
    P3=zeros(1,orderp+1);
    
    patchpos1=post;
 
    patchpos2=post+posdist;
    P3(:)=polyfit(x(1,patchpos1:patchpos2),Wnew(i,patchpos1:patchpos2),orderp);
    Pdd=polyder(P3);
    F(i,j)=(polyval(Pdd,x(1,patchpos2))-polyval(Pdd,x(1,patchpos1)));
    holdpatchpos(:,j)=[patchpos1;patchpos2];
    
    PP(i,:,j)=P3;
    
    subplot(size(Wnew,1),1,i);plot(x(1,patchpos1:patchpos2),polyval(PP(i,:,j),x(1,patchpos1:patchpos2))) %check if polynomials ok
    
    end
    PForce(i,1)=W(i,size(x,2));                %preumond
    PForcen(i,1)=Wnew(i,size(x,2));            %natsiavas
end


T=F';

%Q=Qb*F;
Q=F;
A=[zeros(nummodes,nummodes)  , eye(nummodes) ; -M\K , -M\Z];
B=[zeros(nummodes,numofpatches) ; M\Q];    %for state space control
Bdistba.w=[zeros(nummodes,1) ; M\PForcen]; 
C=eye(2*nummodes);
D=zeros(size(C,1),size(B,2)+1);
Mhold=M; %#ok<NASGU>

%the uncontrollable ss system 4 inputs V and outputs for every mode
sys1=ss(A,[B,Bdistba.w],C,D);
%controllability
Pcontrollability=ctrb(A,B);
unco1=length(A)-rank(Pcontrollability);


%free system - no controllers involved
sys1.OutputName='tau';figure;step(sys1(1:nummodes,numofpatches+1)*force);title('analytic model')
%simulation 
[y,~]=step(sys1(1:nummodes,numofpatches+1)*force);
figure(100);subplot(1,2,1);plot(x,y(size(y,1),1:nummodes)*Wnew);title('analytic modeling')  %natsiavas
ylabel('w [m]');xlabel('L [m]')
%%
%FEMs
%consists of #n FEMs



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
%contraints determination
%eigenvalues & eigenvectors

Kr=Ktot(3:size(Ktot,1),3:size(Ktot,2));Mr=Mtot(3:size(Ktot,1),3:size(Ktot,2));
[vec,valsqr]=eig(Kr,Mr);Eigenvalues=sqrt(diag(valsqr));

%%
%modes
bhma=1/100;                            %points step
le=0:bhma:ltot/n;                      %x-axis position
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
    if i~=size(WH,1); set(gca,'XTick',[]); end ;
end

%%
%Force at the edge again

temp=2;numofpatches=0;
for i=1:floor(n/2)
temp=temp+4;numofpatches=numofpatches+1;
end

b=zeros(2*n,numofpatches);             %determine bs the refering to control actions on evr dof - patches implementation
bw=b(:,1);bw(size(bw,1)-1,1)=1;        %determine the vibration - force at the right edge

temp=2;temp2=1;
for i=1:floor(n/2)
b(temp,temp2)=1;b(temp+2,temp2)=-1;
temp=temp+4;temp2=temp2+1;
end


SQ=sqrt(valsqr);Z=zeros(size(SQ,1),size(SQ,2));
for i=1:2*n
if i==1;Z(i,i)=SQ(i,i)*2*zeta(1,1);end
if i~=1;Z(i,i)=SQ(i,i)*2*zeta(2,1);end
end
Cr=Mr*vec*Z*vec'*Mr;
Ass=[zeros(2*n,2*n),eye(2*n);-inv(Mr)*Kr,-inv(Mr)*Cr];
Bss=[zeros(2*n,2*n);inv(Mr)]*b;
Css=eye(4*n);Dss=0;
Bw=[zeros(2*n,2*n);inv(Mr)]*bw;


sys=ss(Ass,Bw,Css,Dss); 
%figure;step(sys(1:length(Ass)/2)*force);title('FEM model')
clear y t
[y,~]=step(sys(1:length(Ass))*force);

%simulation - SOS for FEMS

yinuse=[zeros(size(y,1),2),y(:,1:2*n)];  %adding 2 more columns zeros in order to overcome size problems
ylast=yinuse(size(yinuse,1),:) ;         %holding the last value - final time
Amplitute=zeros(n,length(le));
temp=0;
for j=1:n                             %for every part

  for i=1:length(le)                   %for part length(plot)
  
   phi1=1-3*(le(i)/(ltot/n))^2+2*(le(i)/(ltot/n))^3;
   phi2=(ltot/n)*((le(i)/(ltot/n))-2*(le(i)/(ltot/n))^2+(le(i)/(ltot/n))^3);
   phi3=3*(le(i)/(ltot/n))^2-2*(le(i)/(ltot/n))^3;
   phi4=(ltot/n)*(-(le(i)/(ltot/n))^2+(le(i)/(ltot/n))^3);
   Amplitute(j,i)=phi1*ylast(1,1+temp)+phi2*ylast(1,2+temp)+phi3*ylast(1,3+temp)+phi4*ylast(1,4+temp);
  
  end 
temp=temp+2;
end


figure(100);subplot(1,2,2);lehold=0;
for i=1:n
plot(le+lehold,Amplitute(i,:))
hold on;
lehold=le(length(le))+lehold;
end
title('FEM analysis');xlabel('L [m]')
%%
%same ? 
Mhold=eye(length(valsqr));
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
AA=[zeros(length(valsqr),length(valsqr))  , eye(length(valsqr)) ; -Mhold\K , -Mhold\Z];
BB=[zeros(length(valsqr),numofpatches) ; Mhold\(vec'*b)];
BBw=[zeros(length(valsqr),1) ; Mhold\(vec'*bw)];
CC=eye(2*length(valsqr));DD=0;
sys3=ss(AA,BBw,CC,DD);
[y,t]=step(sys3*force);
zyy=zeros(size(y,2),1);zy=zeros(size(y,1),length(valsqr));
for i=1:size(y,1)
   zyy=vec*y(i,1:length(valsqr))'; 
   zy(i,:)=zyy';
end
figure;plot(t,zy);title('2 check the preumond method #2 - from fem to disc')
%for sys3 and 
ttemp=t(length(t));
figure(101)
plot(lee',WH'*y(length(t),1:size(y,2)/2)');
%3 way confirmation



%reduced order 

figure(102)
plot(lee',WH(1:n,:)'*y(length(t),1:n)');

figure(103)
plot(lee',WH(n+1:2*n,:)'*y(length(t),n+1:2*n)');


figure(104);temp=1;
for i=1:10
subplot(10,2,temp);plot(t,y(:,i))
if i~=10; set(gca,'XTick',[]); end ;
subplot(10,2,temp+1);plot(t,y(:,i+10))
temp=temp+2;
if i~=10; set(gca,'XTick',[]); end ;
end

