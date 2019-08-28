clear; 

pi=atan(1)*4; ic=sqrt(-1);
damping=10^-7;
Nit=40; %Nit is the maximun iteration number to reach half infinity when solving surface Green's function
L=20;W_pn=10;
Efield=0.1;
delta = 0.2;  deltaterm=diag([delta,delta,-delta,-delta]);
I1=1:4;
I2=5:8;
NN4=4;


t1 = -1.22;
t2 = 3.665;
t3 = -0.205;
t4 = -0.105;
t5 = -0.055;


kx=-pi/3;
ek=exp(ic*kx);
ek1=ek';


nE=0;
cond=0;
for E=-6:0.1:8
    nE=nE+1;
z=E+ic*damping;

H0=0;
H0=[0,t1+t3*ek1,t4*(1+ek1),t5+t2*ek1; t1+t3*ek,0,t2+t5*ek1,t4*(1+ek1);t4*(1+ek),t2+t5*ek,0,t1+t3*ek1; t5+t2*ek,t4*(1+ek),t1+t3*ek,0]+ deltaterm; 
V1=0;
%from the cell y-1 to the cell y
V1=[0,t3*ek1+t1,t4*(1+ek1),0;0,0,0,0;0 0 0 0; 0,t4+t4*ek,t1+t3*ek,0];

V2=0;
 %from the cell y-1 to the cell y
V2=[0 0 0 0;t1+t3*ek,0,0,t4+t4*ek1;t4+t4*ek,0,0,t1+t3*ek1;0 0 0 0;];



 H00=[H0,V1;V2,H0];
 
G0Lead=(z*eye(8)-H00)^-1; 
G011=G0Lead(I1,I1);
G01N=G0Lead(I1,I2);
G0N1=G0Lead(I2,I1);
G0NN=G0Lead(I2,I2);

%%%surface Green's function of infinite electrode%%
ss0=10000;ss=9000;
nit=0;
while (abs(ss0-ss)>1d-9)&(nit<Nit) 
        ss0=ss; nit=nit+1;
G11=G011+G01N*V1*(eye(NN4)-G011*V2*G0NN*V1)^-1*G011*V2*G0N1;  
G2N2N=G0NN+G0N1*V2*(eye(NN4)-G0NN*V1*G011*V2)^-1*G0NN*V1*G01N ; ;

G12N=G01N*V1*(eye(NN4)-G011*V2*G0NN*V1)^-1*G01N;
G2N1=G0N1*V2*(eye(NN4)-G0NN*V1*G011*V2)^-1*G0N1;


G011=G11; G0NN=G2N2N; G01N=G12N; G0N1=G2N1;

ss=sum(sum(abs(G11),1));

end;
nit;
'lead iteration finished';


for i=1:L
H00((i-1)*4+1:i*4,(i-1)*4+1:i*4)=H0;
H00((i-1)*4+1:i*4,(i-1)*4+1:i*4)=H0;
if (i+1<=L)H00(i*4+1:(i+1)*4,(i-1)*4+1:i*4)=V2;end;
if (i+1<=L)H00((i-1)*4+1:i*4,i*4+1:(i+1)*4)=V1;end;
end

for i=1:L
  if (((L - W_pn) / 2 < i )&(i< (L + W_pn) / 2)) 
      ss=(i-(L - W_pn) / 2.0);
      H00((i-1)*4+1:i*4,(i-1)*4+1:i*4)=deltaterm+[ss*Efield,t1+t3*ek1,t4*(1+ek1),t5+t2*ek1;t1+t3*ek,(ss+1d0/2)*Efield,t2+t5*ek1,t4*(1+ek1);t4*(1+ek),t2+t5*ek,(ss+1d0/2)*Efield,t1+t3*ek1;t5+t2*ek,t4*(1+ek),t1+t3*ek,ss*Efield];
  end;
  if i<=(L - W_pn) / 2
       H00((i-1)*4+1:i*4,(i-1)*4+1:i*4)=H0+deltaterm;
  end;
  if i>=(L + W_pn) / 2
       H00((i-1)*4+1:i*4,(i-1)*4+1:i*4)=H0+W_pn*Efield*eye(4)+deltaterm; 
  end;
end;
 G00=(z*eye(L*4)-H00)^-1; 
 
IL=4*(L-1)+1:4*L;
G0=zeros(8,8);
G0(I1,I1)=G00(I1,I1);
G0(I2,I2)=G00(IL,IL);
G0(I1,I2)=G00(I1,IL);
G0(I2,I1)=G00(IL,I1);

selflead=[V2*G0NN*V1,zeros(NN4,NN4);zeros(NN4,NN4),V1*G011*V2];
G0w=selflead*G0;
G0R=G0*(eye(8)-G0w)^-1;



G0A=G0R';
g10=V1*G011*V2;g1=ic*(g10-g10');
g20=V2*G0NN*V1;g2=ic*(g20-g20');


tt=trace(g2*G0R(I1,I2)*g1*G0A(I2,I1));
cond(nE,1)=E;
cond(nE,2)=real(tt);
end
figure;
plot(cond(:,1),cond(:,2),'*-')