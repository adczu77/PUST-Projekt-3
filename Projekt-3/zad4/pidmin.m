function E = pidmin(vector)
Tp=0.5;

Umax = 1;
Umin = -1;
Upp=0;
Ypp=0;
E = 0;

Kr = vector(1);  % 0.6; 0.5; 0.7; 0.8; 0.75;
Ti = vector(2);   % 12; 10; 8; 11; 10;             
Td = vector(3);  % 2.5; 2; 5; 3; 4;    

r2=(Kr*Td)/Tp;
r1=Kr*(Tp/(2*Ti)-2*(Td/Tp)-1);
r0=Kr*(Tp/(2*Ti)+Td/Tp+1);

kk=600;
u(1:kk)=Upp;
y(1:kk)=Ypp;
yzad(1:kk)=-0.3;
yzad(150:kk)=10;
yzad(300:kk)=0.5;
yzad(450:kk)=1.5;
e(1:6)=0;

for k=7:kk
    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));
    e(k)=yzad(k)-y(k); 
    u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
    u(k) = min(u(k), Umax);
    u(k) = max(u(k), Umin);
    E = E + (yzad(k)- y(k))^2;
end