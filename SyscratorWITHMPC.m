function [ sys , Kmpc] = SyscratorWITHMPC( Np , Nc , R , Q , Am , Bm , Cm, SampleTime)
 

A=[Am zeros(size(Am,1),size(Cm,1));Cm*Am eye(size(Cm,1))];
B=[Bm ; Cm*Bm];C=[zeros(size(Cm,1),size(Am,1)) eye(size(Cm,1)) ];
F=zeros(size(Cm,1)*Np,size(Am,1)+size(Cm,1));
phi=zeros(size(Cm,1)*Np,size(Bm,2)*Nc);

for i=1:Np
   temp=i-1;
   F(temp*size(Cm,1)+1:temp*size(Cm,1)+size(Cm,1),:)=C*A^i;
end
for i=1:Np
   temp1=i-1;
   for j=1:Nc
   temp2=j-1;
   exponent=i-j;
   if exponent>=0
       phi(temp1*size(Cm,1)+1:temp1*size(Cm,1)+size(Cm,1),temp2*size(Bm,2)+1:temp2*size(Bm,2)+size(Bm,2))=C*(A^exponent)*B;
   end
   end
end

holdthenextvalue=zeros(size(Bm,2),size(Bm,2)*Nc);holdthenextvalue(1:size(Bm,2),1:size(Bm,2))=eye(size(Bm,2));
Kmpc=(phi'*Q*phi+R)\phi'*F;
Kmpc=holdthenextvalue*Kmpc;

SP=zeros(size(Cm,1)*Np,size(Cm,1));
for i=1:Np
temp1=i-1;
SP(temp1*size(Cm,1)+1:temp1*size(Cm,1)+size(Cm,1),:)=eye(size(Cm,1));
end;

%taus sp =1 and time derivative of each tau =0 respectively 
%Ysp=ones(size(Cm,1),1);Ysp((size(Cm,1)/2)+1:size(Cm,1))=zeros(size(Cm,1)/2,1);

%set points without disturbances
%Ky=(phi'*Q*phi+R)\phi'*SP*Ysp;
Ky=(phi'*Q*phi+R)\phi'*SP;
Ky=holdthenextvalue*Ky;



% Es=[zeros(16,8);eye(8);zeros(8,8)];
% sysplant=ss(A,[Es,B,zeros(size(Kmpc'))],eye(size(A,1)),0,SampleTime);
% sysmodel=ss(A,[zeros(size(Es)),B,Kmpc'],eye(size(A,1)),0,SampleTime);
sys=ss(A-B*Kmpc,B*Ky,eye(length(A)),0,SampleTime);

end

