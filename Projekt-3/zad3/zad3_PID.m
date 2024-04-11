% Inicjalizacja zmiennych:
clear;
Tp=0.5;

Umax = 1;
Umin = -1;
Upp=0;
Ypp=0;
E = 0;

% Ustawienie parametrów:
Kr = 0.1;  
Ti = 4;             
Td = 2;    

% Wyliczenie zmiennych r:
r2=(Kr*Td)/Tp;
r1=Kr*(Tp/(2*Ti)-2*(Td/Tp)-1);
r0=Kr*(Tp/(2*Ti)+Td/Tp+1);

% Ustalenie warunków początkowych
kk=600;
u(1:kk)=Upp;
y(1:kk)=Ypp;
yzad(1:kk)=-0.3;
yzad(150:kk)=10;
yzad(300:kk)=0.5;
yzad(450:kk)=1.5;
e(1:6)=0;

% Algorytm PID:
for k=7:kk
    % Pomiar wyjścia:
    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));
    % Wyliczenie uchybu:
    e(k)=yzad(k)-y(k); 
    % Wyliczenie sterowania:
    u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
    % Uwzględnienie ograniczeń:
    u(k) = min(u(k), Umax);
    u(k) = max(u(k), Umin);
    % Wyliczanie błędu:
    E = E + (yzad(k)- y(k))^2;
end

% Pokazanie przebiegu:
figure; 
stairs((1:kk),u); 
title("u, E="+ num2str(E)); 
xlabel('k'); 
ylabel('u');
%print('pust5_ufmincon.png','-dpng','-r400');

figure;
hold on;
stairs(1:kk,y);  
stairs(1:kk,yzad,':'); 
title("yzad, y, E="+ num2str(E)); 
xlabel('k'); 
ylabel('y');
legend('y','yzad');
%print('pust5_yfmincon.png','-dpng','-r400');