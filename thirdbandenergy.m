function energy3=thirdbandenergy(kly)
global kx delta vpn Efermi
i=sqrt(-1);
t1 = -1.22;
t2 = 3.665;
t3 = -0.205;
t4 = -0.105;
t5 = -0.055;


klx=kx;
%for klx=0:2*pi/30:2*pi;
band=0;


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
H=H+vpn*eye(4);
band=sort(eig(H));
energy3=abs(band(3)-Efermi);