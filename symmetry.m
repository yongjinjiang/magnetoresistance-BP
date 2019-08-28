
delta = 0;  
delta1=0.3;
delta2=0.6;
delta3=0.6;
delta4=0.3;

ic=sqrt(-1);
t1 = -1.22;
t2 = 3.665;
t3 = -0.205;
t4 = -0.105;
t5 = -0.055;


kx0=-pi/3;ky0=pi/4;
kx=kx0;ky=ky0;
ek=exp(ic*kx); ey=exp(ic*ky);ey1=ey';
ek1=ek';
H0=[delta1,(t1+t3*ek1)*(1+ey1),t4*(1+ek1)*(1+ey1),t5+t2*ek1; 
   0,delta2,t2+t5*ek1,t4*(1+ek1)*(1+ey);
   0,0,delta3,(t1+t3*ek1)*(1+ey); 
    0,0,0,delta4]; 
H0=H0+H0';

P=[0,1,0,0;
    1,0,0,0;
    0,0,0,ek;
    0,0,ek,0];
P1=P';
M=[ey1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,ey1];
M1=M';

PP=[0,0,0,ek;
    0,0,ek*ey,0;
    0,ek*ey,0,0;
    ek,0,0,0
     ];
 PP1=PP';

kx=-kx0;ky=-ky0;
ek=exp(ic*kx); ey=exp(ic*ky);ey1=ey';
ek1=ek';
H1=[delta1,(t1+t3*ek1)*(1+ey1),t4*(1+ek1)*(1+ey1),t5+t2*ek1; 
   0,delta2,t2+t5*ek1,t4*(1+ek1)*(1+ey);
   0,0,delta3,(t1+t3*ek1)*(1+ey); 
    0,0,0,delta4]; 
H1=H1+H1';

P1symmetry=P1*H0*P-H1;
PP1symmetry=PP1*H0*PP-H1;

kx=kx0;ky=-ky0;
ek=exp(ic*kx); ey=exp(ic*ky);ey1=ey';
ek1=ek';
H2=[delta1,(t1+t3*ek1)*(1+ey1),t4*(1+ek1)*(1+ey1),t5+t2*ek1; 
   0,delta2,t2+t5*ek1,t4*(1+ek1)*(1+ey);
   0,0,delta3,(t1+t3*ek1)*(1+ey); 
    0,0,0,delta4]; 
H2=H2+H2';

Msymmetry=M1*H0*M-H2