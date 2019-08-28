clear; 
damping=10^-5;
pi=atan(1)*4; ic=sqrt(-1);


N=20;      %N is the width of the ribbon.
NN4=N*1*2; %number of states of each layer
NN8=NN4*2;  % sites needed for each iteration for Greens' function of the disordered system. 
             %number 4 and 8 are due to previous project where each site has four internal states 

I1=1:NN4; I2=NN4+1:NN8;

NL=198; %NL*2 is the length of the disordered system. Later we periodically add these systems into infinite(not half-infinite!) ribbon 
Nit=40; %Nit is the maximun iteration number to reach half infinity when solving surface Green's function

W=0; U=0;
M=-10;
E=14;
z=E+ic*damping;

N_disorder=50; N_W=30;

a=5d0; 
A=364.5;
B=-686;
C=0;
D=-512;

Es=C+M-4*(B+D)/a^2;
Ep=C-M-4*(D-B)/a^2;
Vss=(B+D)/a^2;
Vpp=(D-B)/a^2;
Vsp=-ic*A/(2*a);



sigma_x=[0,1;1,0];
sigma_y=[0,-ic;ic,0];
sigma_z=[1,0;0,-1];



i=0;
for nx=1:N   
    for a=1:2
        i=i+1;
        ind(a,nx)=i;
    end;
end



H0=0;
for nx=1:N
         %%%%hoping at x direction    open boundary 
           
          if (nx+1<=N)
           p=mod(nx+1,N+1)+(1-(-1)^floor((nx+1)/(N+1)))/2;
           H0(ind(1:2,nx),ind(1:2,p))=-Vsp*sigma_x;
           H0(ind(1,nx),ind(1,p))=Vss; 
           H0(ind(2,nx),ind(2,p))=Vpp;
           H0(ind(1:2,nx)+NN4,ind(1:2,p)+NN4)=-Vsp*sigma_x;
           H0(ind(1,nx)+NN4,ind(1,p)+NN4)=Vss; 
           H0(ind(2,nx)+NN4,ind(2,p)+NN4)=Vpp;
           end 
           if (nx-1>=1)   
           p=nx-1+N*(1-(-1)^floor(((nx-1)-0.1)/N))/2;
           H0(ind(1:2,nx),ind(1:2,p))=Vsp*sigma_x;
           H0(ind(1,nx),ind(1,p))=Vss; 
           H0(ind(2,nx),ind(2,p))=Vpp;
           H0(ind(1:2,nx)+NN4,ind(1:2,p)+NN4)=Vsp*sigma_x;
           H0(ind(1,nx)+NN4,ind(1,p)+NN4)=Vss; 
           H0(ind(2,nx)+NN4,ind(2,p)+NN4)=Vpp;
          end 
 end;


  for nx=1:N
        H0(ind(1:2,nx),ind(1:2,nx))=zeros(2,2);
        H0(ind(1,nx),ind(1,nx))=H0(ind(1,nx),ind(1,nx))+Es-10;
        H0(ind(2,nx),ind(2,nx))=H0(ind(2,nx),ind(2,nx))+Ep-10;
        
        H0(ind(1:2,nx)+NN4,ind(1:2,nx)+NN4)=zeros(2,2);
        H0(ind(1,nx)+NN4,ind(1,nx)+NN4)=H0(ind(1,nx)+NN4,ind(1,nx)+NN4)+Es-10;
        H0(ind(2,nx)+NN4,ind(2,nx)+NN4)=H0(ind(2,nx)+NN4,ind(2,nx)+NN4)+Ep-10;
end
for nx=1:N
 
        H0(ind(1:2,nx),ind(1:2,nx)+NN4)=-Vsp*sigma_y;
         H0(ind(1,nx),ind(1,nx)+NN4)=Vss;
          H0(ind(2,nx),ind(2,nx)+NN4)=Vpp;
         
         H0(ind(1:2,nx)+NN4,ind(1:2,nx))=Vsp*sigma_y;
         H0(ind(1,nx)+NN4,ind(1,nx))=Vss;
          H0(ind(2,nx)+NN4,ind(2,nx))=Vpp; 
end
V1=H0(I1,I2);
V2=V1';

G0Lead=(z*eye(NN8)-H0)^-1; 


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
nit
'lead iteration finished';



for nW=1:N_W
    W=10*(nW-1);
  cond=0;cond(N_disorder)=0;

  for n_disorder=1:N_disorder
H0=zeros(NN8); %two layers of section-sites constitute a iterative unit 
V1=0;V2=0;

for nx=1:N
        H0(ind(1:2,nx),ind(1:2,nx))=W*(rand-0.5)*eye(2)+U*(rand-0.5)*sigma_x;
        H0(ind(1,nx),ind(1,nx))=H0(ind(1,nx),ind(1,nx))+Es;
        H0(ind(2,nx),ind(2,nx))=H0(ind(2,nx),ind(2,nx))+Ep;
        
        H0(ind(1:2,nx)+NN4,ind(1:2,nx)+NN4)=W*(rand-0.5)*eye(2)+U*(rand-0.5)*sigma_x;
        H0(ind(1,nx)+NN4,ind(1,nx)+NN4)=H0(ind(1,nx)+NN4,ind(1,nx)+NN4)+Es;
        H0(ind(2,nx)+NN4,ind(2,nx)+NN4)=H0(ind(2,nx)+NN4,ind(2,nx)+NN4)+Ep;
end



for nx=1:N
         %%%%hoping at x direction    open boundary 
           
          if (nx+1<=N)
           p=mod(nx+1,N+1)+(1-(-1)^floor((nx+1)/(N+1)))/2;
           H0(ind(1:2,nx),ind(1:2,p))=-Vsp*sigma_x;
           H0(ind(1,nx),ind(1,p))=Vss; 
           H0(ind(2,nx),ind(2,p))=Vpp;
           H0(ind(1:2,nx)+NN4,ind(1:2,p)+NN4)=-Vsp*sigma_x;
           H0(ind(1,nx)+NN4,ind(1,p)+NN4)=Vss; 
           H0(ind(2,nx)+NN4,ind(2,p)+NN4)=Vpp;
           end 
           if (nx-1>=1)   
           p=nx-1+N*(1-(-1)^floor(((nx-1)-0.1)/N))/2;
           H0(ind(1:2,nx),ind(1:2,p))=Vsp*sigma_x;
           H0(ind(1,nx),ind(1,p))=Vss; 
           H0(ind(2,nx),ind(2,p))=Vpp;
           H0(ind(1:2,nx)+NN4,ind(1:2,p)+NN4)=Vsp*sigma_x;
           H0(ind(1,nx)+NN4,ind(1,p)+NN4)=Vss; 
           H0(ind(2,nx)+NN4,ind(2,p)+NN4)=Vpp;
          end 
 end;


H00=H0; %for later interaion use



for nx=1:N
 
        H0(ind(1:2,nx),ind(1:2,nx)+NN4)=-Vsp*sigma_y;
         H0(ind(1,nx),ind(1,nx)+NN4)=Vss;
          H0(ind(2,nx),ind(2,nx)+NN4)=Vpp;
         
         H0(ind(1:2,nx)+NN4,ind(1:2,nx))=Vsp*sigma_y;
         H0(ind(1,nx)+NN4,ind(1,nx))=Vss;
          H0(ind(2,nx)+NN4,ind(2,nx))=Vpp; 
end
%V1=H0(ind(1,1,1:4),ind(1,1,1:4)+4NN);V2=V1';
V1=H0(I1,I2);
V2=V1';
G0=(z*eye(NN8)-H0)^-1; 


%%%%%%%%%%%calculate the surface Greens' function of the disordered slab iteratively%%%%%%%%%%%%%% 

for i=1:NL
    H0=H00;%update diagonal part of two-plane slab
  for nx=1:N
        H0(ind(1:2,nx),ind(1:2,nx))=W*(rand-0.5)*eye(2)+U*(rand-0.5)*sigma_x;
        H0(ind(1,nx),ind(1,nx))=H0(ind(1,nx),ind(1,nx))+Es;
        H0(ind(2,nx),ind(2,nx))=H0(ind(2,nx),ind(2,nx))+Ep;
        
        H0(ind(1:2,nx)+NN4,ind(1:2,nx)+NN4)=W*(rand-0.5)*eye(2)+U*(rand-0.5)*sigma_x;
        H0(ind(1,nx)+NN4,ind(1,nx)+NN4)=H0(ind(1,nx)+NN4,ind(1,nx)+NN4)+Es;
        H0(ind(2,nx)+NN4,ind(2,nx)+NN4)=H0(ind(2,nx)+NN4,ind(2,nx)+NN4)+Ep;
end


    %% glue the two planes by Green's function of old slab  
     H0(I1,I1)=H0(I1,I1)+V1*G0(I1,I1)*V2; 
     H0(I2,I2)=H0(I2,I2)+V2*G0(I2,I2)*V1; 
     H0(I1,I2)= V1*G0(I1,I2)*V1;
     H0(I2,I1)= V2*G0(I2,I1)*V2;
     G0=(z*eye(NN8)-H0)^-1; 
end;




%%%%%%%%%%%calculate Gpq%%%%%%%%%%%%%%
%H00=H00;

selflead=[V2*G0NN*V1,zeros(NN4,NN4);zeros(NN4,NN4),V1*G011*V2];
G0w=selflead*G0;
G0R=G0*(eye(NN8)-G0w)^-1;
G0A=G0R';
g10=V1*G011*V2;g1=ic*(g10-g10');
g20=V2*G0NN*V1;g2=ic*(g20-g20');


tt=trace(g2*G0R(I1,I2)*g1*G0A(I2,I1));
nW
n_disorder
cond(n_disorder)=real(tt);
sum(cond)/n_disorder
end;

conduct(nW,1)=W;
conduct(nW,2)=sum(cond)/N_disorder;
   dv=0;for i=1:N_disorder; dv=dv+(cond(i)-conduct(nW,2))^2;end;dv=dv/N_disorder; dv=sqrt(dv);
conduct(nW,3)=dv;
end;
conduct2=conduct;
figure;
plot(conduct(:,1),conduct(:,2),'g*',conduct(:,1),conduct(:,3),'r+')
