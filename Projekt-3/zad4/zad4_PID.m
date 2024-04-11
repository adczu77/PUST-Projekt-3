clear;
Tp=0.5;

Umax = 1;
Umin = -1;
Upp=0;
Ypp=0;
E = 0;

Kr = 0.1471;  % 0.1; 0.15; 0.1; 0.1; 0.05;
Ti = 2.5621;   % 6; 6; 8; 4; 4;             
Td = 0.0972;  % 1; 1.5; 0.5; 3; 2;    

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
e(1:kk)=0;

for k=7:kk
    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));
    e(k)=yzad(k)-y(k); 
    u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
    u(k) = min(u(k), Umax);
    u(k) = max(u(k), Umin);
    E = E + (yzad(k)- y(k))^2;
end


figure; 
stairs((1:kk),u); 
title("u, E="+ num2str(E)); 
xlabel('k'); 
ylabel('u');
print('pust_pid_u_fmincon.png','-dpng','-r400');

figure;
hold on;
stairs(1:kk,y);  
stairs(1:kk,yzad,':'); 
title("yzad, y, E="+ num2str(E)); 
xlabel('k'); 
ylabel('y');
legend('y','yzad');
print('pust_pid_y_fmincon.png','-dpng','-r400');