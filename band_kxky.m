

t1 = -1.22;
t2 = 3.665;
t3 = -0.205;
t4 = -0.105;
t5 = -0.055;
delta = 0.5;


klx=-pi/2;
band=0;
nkx=0;
for klx=0:2*pi/100:2*pi- 2*pi/100
   nkx=nxk+1;
   nky=0;
for kly=0:2*pi/100:2*pi- 2*pi/100
    H=0;nky=nky+1;
H12 = t1*(1+exp(-i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(-i*kly));  
H13 = t4*(1+exp(-i*kly)+exp(-i*klx)+exp(-i*klx)*exp(-i*kly));         
H14 = t5 + t2*exp(-i*klx);                                            
H23 = t2 + t5*exp(-i*klx);                                           
H24 = t4*(1+exp(i*kly)+exp(-i*klx)+exp(-i*klx)*exp(i*kly));           
H34 = t1*(1+exp(i*kly)) + t3*(exp(-i*klx)+exp(-i*klx)*exp(i*kly));     
H = [0 H12 H13 H14;
    conj(H12) 0 H23 H24;
    conj(H13) conj(H23) delta H34;
    conj(H14) conj(H24) conj(H34) delta];
band(nkx,nky,1:2)=[klx,kly];
band(nkx,nky,2:5)=eig(H)';
end;
band(:,:,3)