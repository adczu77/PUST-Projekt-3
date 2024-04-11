clear;
Tp=0.5;

Umax = 1;
Umin = -1;
Upp=0;
Ypp=0;
E = 0;
liczba_regulatorow=5;
szerokosc=Umax-Umin;
if liczba_regulatorow==2
    Kr(1:liczba_regulatorow)=[4 0.05];
    Ti(1:liczba_regulatorow) = [2 4];              
    Td(1:liczba_regulatorow) = [0.1 0.4];
    bok = szerokosc/(liczba_regulatorow+1);
    mf = [fismf("trapmf",[-100 -1 -0.8 -0.6]);
      fismf("trapmf",[-0.8 -0.6 1 100]);
      ];
end
if liczba_regulatorow==3
    Kr(1:liczba_regulatorow)=[3 0.15 0.035];
    Ti(1:liczba_regulatorow) = [2 1.75 4];              
    Td(1:liczba_regulatorow) = [0.5 0.9 0.4];
    bok = szerokosc/(liczba_regulatorow+1);
    mf = [fismf("trapmf",[-100 -1 -0.8 0]);
      fismf("trapmf",[-0.8 0 0 0.2]);
      fismf("trapmf",[0 0.2 1 100]);
];
end
if liczba_regulatorow==4
    Kr(1:liczba_regulatorow)=[2 0.75 0.1 0.03];
    Ti(1:liczba_regulatorow) = [1.2 1.95 3.2 4];              
    Td(1:liczba_regulatorow) = [0.6 0.8 0.5 0.7];
    bok = szerokosc/(liczba_regulatorow+1);
    mf = [fismf("trapmf",[-100 -1 -0.6 -0.4]);
      fismf("trapmf",[-0.6 -1+bok*2 -1+bok*2 0]);
      fismf("trapmf",[-0.2 -1+bok*3 -1+bok*3 0.3]);
      fismf("trapmf",[0.2 0.3 1 100]);
];
end
if liczba_regulatorow==5
    Kr(1:liczba_regulatorow)=[3 1.5 0.15 0.05 0.01];
    Ti(1:liczba_regulatorow) = [2 1.9 1.7 2.5 2.7];              
    Td(1:liczba_regulatorow) = [0.3 0.4 0.8 0.9 1.5];
    bok = szerokosc/(liczba_regulatorow+1);
    mf = [fismf("trapmf",[-100 -1 -0.8 -0.6]);
      fismf("trapmf",[-0.8 -1+bok*2 -1+bok*2 0]);
      fismf("trapmf",[-0.333 -1+bok*3 -1+bok*3 0.2]);
      fismf("trapmf",[0 -1+bok*4 -1+bok*4 0.4]);
      fismf("trapmf",[0.333 0.4 1 100]);
      ];
end
x = -1:0.001:1;
mi = evalmf(mf,x);

figure;
plot(x,mi)
xlabel('Input (x)')
ylabel('Membership value (y)')

r0(1:liczba_regulatorow)=0;
r1(1:liczba_regulatorow)=0;
r2(1:liczba_regulatorow)=0;

for i=1:liczba_regulatorow
    r0(i)= Kr(i)*(1+ Tp/(2*Ti(i))+(Td(i)/Tp));
    r1(i)= Kr(i)*(Tp/(2*Ti(i))-(2*Td(i)/Tp)-1);
    r2(i)= Kr(i)*Td(i)/Tp;
end

kk=600;
u(1:kk)=Upp;
y(1:kk)=Ypp;
yzad(1:kk)=-0.3;
yzad(150:kk)=10;
yzad(300:kk)=0.5;
yzad(450:kk)=1.5;
e(1:6)=0;
du(1:liczba_regulatorow)=0;

for k=7:kk
    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));
    e(k)=yzad(k)-y(k);
    mi=evalmf(mf,u(k-1));
    for i =1:liczba_regulatorow
        du(i)=r2(i)*e(k-2)+r1(i)*e(k -1)+r0(i)*e(k);
    end
    if liczba_regulatorow==2
        du_suma=(du(1)*mi(1)+du(2)*mi(2))/(mi(1)+mi(2));
        u(k)=u(k-1)+du_suma;
    end
    if liczba_regulatorow==3
        du_suma=(du(1)*mi(1)+du(2)*mi(2)+du(3)*mi(3))/(mi(1)+mi(2)+mi(3));
        u(k)=u(k-1)+du_suma;
    end
    if liczba_regulatorow==4
        du_suma=(du(1)*mi(1)+du(2)*mi(2)+du(3)*mi(3)+du(4)*mi(4))/(mi(1)+mi(2)+mi(3)+mi(4));
        u(k)=u(k-1)+du_suma;
    end
    if liczba_regulatorow==5
        du_suma=(du(1)*mi(1)+du(2)*mi(2)+du(3)*mi(3)+du(4)*mi(4)+du(5)*mi(5))/(mi(1)+mi(2)+mi(3)+mi(4)+mi(5));
        u(k)=u(k-1)+du_suma;
    end
    u(k) = min(u(k), Umax);
    u(k) = max(u(k), Umin);
    E = E + (yzad(k)- y(k))^2;
end

figure; 
stairs((1:kk),u); 
title("u, E="+ num2str(E)); 
xlabel('k'); 
ylabel('u');
print('pust_pid_u_5_regulatory_po.png','-dpng','-r400');

figure;
hold on;
stairs(1:kk,y);  
stairs(1:kk,yzad,':'); 
title("yzad, y, E="+ num2str(E)); 
xlabel('k'); 
ylabel('y');
legend('y','yzad');
print('pust_pid_y_5_regulatory_po.png','-dpng','-r400');
