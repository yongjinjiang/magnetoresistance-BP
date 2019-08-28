
i=sqrt(-1);
t1 = -1.22;
t2 = 3.665;
t3 = -0.205;
t4 = -0.105;
t5 = -0.055;
delta = 1;

klx=-2*pi*0.1*0;
%for klx=0:2*pi/30:2*pi;
band=0;
Efermi=1;

   nky=0;
for kly=-pi:2*pi/100:pi- 2*pi/100
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
figure;
plot(band(:,1),band(:,2:5));hold on; 
title(num2str(klx/(2*pi)));
