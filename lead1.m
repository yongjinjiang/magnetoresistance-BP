function [kyin,Vin,Vout]=lead1()
global kx delta  Efermi veps
i=sqrt(-1);
t1 = -1.22;
t2 = 3.665;
t3 = -0.205;
t4 = -0.105;
t5 = -0.055;

klx=kx;
%for klx=0:2*pi/30:2*pi;
band=0;


   nky=0;
for kly=-pi:2*pi/500:pi- 2*pi/500
    H=0;nky=nky+1;
H12 = t1*(1+exp(-i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(-i*kly));  
H13 = t4*(1+exp(-i*kly)+exp(-i*klx)+exp(-i*klx)*exp(-i*kly));         
H14 = t5 + t2*exp(-i*klx);                                            
H23 = t2 + t5*exp(-i*klx);                                           
H24 = t4*(1+exp(i*kly)+exp(-i*klx)+exp(-i*klx)*exp(i*kly));           
H34 = t1*(1+exp(i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(i*kly));     
H = [delta H12 H13 H14;
    conj(H12) delta H23 H24;
    conj(H13) conj(H23) -delta H34;
    conj(H14) conj(H24) conj(H34) -delta];
band(nky,1)=kly;
band(nky,2:5)=sort(eig(H))';
end;
%figure;
%plot(band(:,1),band(:,2:5));hold on; 
%title(num2str(klx/(2*pi)));
bb=min(band(:,4));

aa=Efermi-band(:,4);[bb,Ib]=sort(abs(aa));

[kyin0,in]=max(band(Ib(1:2),1));  % pick the forward propagating mode

options=optimset;
options.TolFun=10^-12;
[kyin,fval] = fminsearch('thirdbandenergy1',kyin0,options);
%Efermi= band(Ib(in),4);

ly=1;i=sqrt(-1);kly=kyin;
H12 = i*(t1*(-ly*exp(-i*kly)) + t3*(-ly*exp(-i*klx).*exp(-i*kly)))   ;
H13 = i*t4*(-ly*exp(-i*kly)-ly*exp(-i*klx).*exp(-i*kly))            ; 
H14 = 0                                                           ;
H23 = 0                                                             ;
H24 = i*t4*(ly*exp(i*kly)+ly*exp(-i*klx).*exp(i*kly))               ;
H34 = i*(t1*(ly*exp(i*kly)) + t3*(ly*exp(-i*klx).*exp(i*kly)))      ;


Hdky(1,1) = 0;
Hdky(1,2) = H12;
Hdky(1,3) = H13;
Hdky(1,4) = H14;
Hdky(2,1) = conj(H12);
Hdky(2,2)= 0;
Hdky(2,3)= H23;
Hdky(2,4) = H24;
Hdky(3,1) = conj(H13);
Hdky(3,2) = conj(H23);
Hdky(3,3) = 0;
Hdky(3,4)= H34;
Hdky(4,1)= conj(H14);
Hdky(4,2) = conj(H24);
Hdky(4,3)= conj(H34);
Hdky(4,4) = 0;



kly=kyin+veps;

H=0;
H12 = t1*(1+exp(-i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(-i*kly));  
H13 = t4*(1+exp(-i*kly)+exp(-i*klx)+exp(-i*klx)*exp(-i*kly));         
H14 = t5 + t2*exp(-i*klx);                                            
H23 = t2 + t5*exp(-i*klx);                                           
H24 = t4*(1+exp(i*kly)+exp(-i*klx)+exp(-i*klx)*exp(i*kly));           
H34 = t1*(1+exp(i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(i*kly));     
H = [delta H12 H13 H14;
    conj(H12) delta H23 H24;
    conj(H13) conj(H23) -delta H34;
    conj(H14) conj(H24) conj(H34) -delta];
newE=sort(eig(H));
velocity10=(newE(3)-Efermi)/veps;

%options=optimset;
%options.TolFun=10^-12;
%[kyin,fval] = fminsearch('thirdbandenergy',kyin0,options);

kly=kyin;
H=0;
H12 = t1*(1+exp(-i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(-i*kly));  
H13 = t4*(1+exp(-i*kly)+exp(-i*klx)+exp(-i*klx)*exp(-i*kly));         
H14 = t5 + t2*exp(-i*klx);                                            
H23 = t2 + t5*exp(-i*klx);                                           
H24 = t4*(1+exp(i*kly)+exp(-i*klx)+exp(-i*klx)*exp(i*kly));           
H34 = t1*(1+exp(i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(i*kly));     
H = [delta H12 H13 H14;
    conj(H12) delta H23 H24;
    conj(H13) conj(H23) -delta H34;
    conj(H14) conj(H24) conj(H34) -delta];

[v,D]=eig(H);D=diag(D);
Vin=v(:,3);
velocity1=real(Vin'*Hdky*Vin);

%'Efermi-D(3)=check===',Efermi-D(3);

kly=-kyin;
H=0;
H12 = t1*(1+exp(-i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(-i*kly));  
H13 = t4*(1+exp(-i*kly)+exp(-i*klx)+exp(-i*klx)*exp(-i*kly));         
H14 = t5 + t2*exp(-i*klx);                                            
H23 = t2 + t5*exp(-i*klx);                                           
H24 = t4*(1+exp(i*kly)+exp(-i*klx)+exp(-i*klx)*exp(i*kly));           
H34 = t1*(1+exp(i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(i*kly));     
H = [delta H12 H13 H14;
    conj(H12) delta H23 H24;
    conj(H13) conj(H23) -delta H34;
    conj(H14) conj(H24) conj(H34) -delta];

[v,D]=eig(H);D=diag(D);
Vout=v(:,3);
%'Efermi-D(3)=check===',Efermi-D(3);

[Vin,Vout];
ey=exp(i*kyin);ey1=ey';
M=[ey1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,ey1];
M1=M';
Vin./(M*Vout);

Vin=Vin/sqrt(velocity1);
Vout=Vout/sqrt(velocity1);

