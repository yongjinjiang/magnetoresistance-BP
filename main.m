%this program only calculate for P1-P2 junction where Fermi energy lies in
%the third level(lower conduction band)
%clear
global kx delta vpn Efermi veps
ic=sqrt(-1);
veps=0.3848*1d-3;
veps=1d-4;
ic=sqrt(-1); t1 = -1.22; t2 = 3.665; t3 = -0.205; t4 = -0.105; t5 = -0.055;
L=1000; W_pn=L/10;
pottype=2;
kx=-pi/16*0;
delta=0.2;


vpn=0.1;



Efermi=0.85;
deltaterm=zeros(4,4); 
deltaterm(1,1)=delta;   deltaterm(2,2)=delta;
deltaterm(3,3)=-delta;  deltaterm(4,4)=-delta;

[kyin,Vin,Vout]=lead1(); %note Efermi has been changed a little bit
[ky2in,V2in,V2out]=lead2(); 

for ny=1:L
    for ni=1:4  
            id(ni,ny)=ni+(ny-1)*4;
    end;
end;
ek=exp(ic*kx);ek1=ek';
H0=[0,t1+t3*ek1,t4*(1+ek1),t5+t2*ek1; t1+t3*ek,0,t2+t5*ek1,t4*(1+ek1);t4*(1+ek),t2+t5*ek,0,t1+t3*ek1; t5+t2*ek,t4*(1+ek),t1+t3*ek,0]+ deltaterm; 
V1=0;
%from the cell y-1 to the cell y  from 2,3 to 1,4
V1=[0,t3*ek1+t1,t4*(1+ek1),0;0,0,0,0;0 0 0 0; 0,t4+t4*ek,t1+t3*ek,0];

V2=0;
 %from the cell y to the cell y-1
V2=[0 0 0 0;t1+t3*ek,0,0,t4+t4*ek1;t4+t4*ek,0,0,t1+t3*ek1;0 0 0 0;];

%%%%%%%H00 is the Hamiltonina for the whole system%%%%%%%%%%%%%%%%
for ny=1:L
   H00(id(:,ny),id(:,ny))=H0;
end;
for ny=1:L
    if (pottype==1)
       if (ny>(L-W_pn)/2)&(ny<(L+W_pn)/2);
         ss=ny-(L-W_pn)/2;
         H00(id(1,ny),id(1,ny))=vpn/W_pn*ss+H00(id(1,ny),id(1,ny));        H00(id(4,ny),id(4,ny))=vpn/W_pn*ss+H00(id(4,ny),id(4,ny));
         H00(id(2,ny),id(2,ny))=vpn/W_pn*(ss+0.5)+H00(id(2,ny),id(2,ny));  H00(id(3,ny),id(3,ny))=vpn/W_pn*(ss+0.5)+H00(id(3,ny),id(3,ny));
       end;
       if (ny>=(L+W_pn)/2)
         H00(id(:,ny),id(:,ny))=H00(id(:,ny),id(:,ny))+vpn*eye(4);
       end;
      else
       ss=ny-(L-W_pn)/2;
      H00(id(1,ny),id(1,ny))=vpn/2*(1+tanh(ss/W_pn))+H00(id(1,ny),id(1,ny));        H00(id(4,ny),id(4,ny))=vpn/2*(1+tanh(ss/W_pn))+H00(id(4,ny),id(4,ny));
      H00(id(2,ny),id(2,ny))=vpn/2*(1+tanh((ss+0.5)/W_pn))+H00(id(2,ny),id(2,ny));  H00(id(3,ny),id(3,ny))=vpn/2*(1+tanh((ss+0.5)/W_pn))+H00(id(3,ny),id(3,ny));
    end;
end;



for ny=1:L-1
H00(id(:,ny+1),id(:,ny))=V1;
H00(id(:,ny),id(:,ny+1))=V2;
end

%for ny=


%%%%%%%%%%add more terms in H00%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%H00(100,100)=100;
H02=H0+vpn*eye(4);
A=H00-Efermi*eye(4*L);  A1=A^-1;

C=zeros(2,2);
v0=0;v0=Efermi*Vout-H0*Vout; C(1,1)= v0(2); % for the downlead, we choose the site 2 
v0=0;v0=Efermi*V2out-H02*V2out; C(2,2)= v0(1); % for the upper lead, we choose site 1 for constructing independent equation

D=zeros(2,4*L);
D(1,id(:,1))=V2(2,:);
D(2,id(:,L))=V1(1,:);

B=zeros(4*L,2);
B(id(:,1),1)= -V1*Vout;  B(id(:,L),2)=-V2*V2out;
%%%%%%%%%%%%%from left to right%%%%%%%%%%%%% 
sL=1;sR=0;

b=zeros(4*L,1);
b(id(:,1))=-V1*sL*Vin;   %left(or down)lead to central region
b(id(:,L))=-V2*sR*V2in;  % right(or upper) lead to central region ,sL and sR means where to cincident


d=zeros(2,1);
v0=0;v0=(H0-Efermi*eye(4))*Vin*sL; d(1)=v0(2); 
v0=0;v0=(H02-Efermi*eye(4))*V2in*sR; d(2)=v0(1); 
 

g=d+D*A1*b;
M=C-D*A1*B;
r1=M^-1*g;

Efermi

abs(r1).^2

%%check:
psai=A1*(b+B*r1);
aa=A*psai;check1=aa(id(:,1))+V1*(Vin*sL+r1(1)*Vout);
aa=A*psai;check2=aa(id(:,L))+V2*(V2in*sR+r1(2)*V2out);


 ek=exp(ic*kx);
 ek1=exp(-ic*kx);
 
 for ny=1:L-1
     wfn=psai(id(:,ny));
     wfn1=psai(id(:,ny+1));
 
current0=0;
current0=t2* (wfn(2))'*wfn(3)+t5* (wfn(1))'*wfn(4)+ t4* (wfn(2))'*wfn(4)+t4* (wfn(1))'*wfn(3)+ (t3* ek*(wfn(2))'*wfn(1)+ t3* ek*(wfn(4))'*wfn(3))...
        +t4*((wfn1(1))'*wfn(3) +(wfn(2))'*wfn1(4)) +t3*ek*((wfn1(4))'*wfn(3) + (wfn(2))'*wfn1(1));
current1=current0;
current2(ny)=(current1-(current1)')/ic;

 end 
figure;plot(1:L-1,current2);
figure;plot(1:L-1,abs(psai(id(1,1:L-1))),'*-')

% sum(current);
% 
% stop



































%%%%%%%%%%%%%from right to left%%%%%%%%%%%%% 
sL=0;sR=1;

b=zeros(4*L,1);
b(id(:,1))=-V1*sL*Vin;   %left(or down)lead to central region
b(id(:,L))=-V2*sR*V2in;  % right(or upper) lead to central region ,sL and sR means where to cincident


d=zeros(2,1);
v0=0;v0=(H0-Efermi*eye(4))*Vin*sL; d(1)=v0(2); 
v0=0;v0=(H02-Efermi*eye(4))*V2in*sR; d(2)=v0(1); 
 

g=d+D*A1*b;
M=C-D*A1*B;
r2=M^-1*g;

%%check:
psai=A1*(b+B*r2);
aa=A*psai;check1=aa(id(:,1))+V1*(Vin*sL+r2(1)*Vout);
aa=A*psai;check2=aa(id(:,L))+V2*(V2in*sR+r2(2)*V2out);
%r2

S=[r1,r2];S'*S;
accuracy=r1'*r1;

