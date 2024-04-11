% Inicjalizacja zmiennych:
clear;
Tp=0.5;

kk=600;
Upp=0;
Ypp=0;
Umax = 1;
Umin = -1;
E = 0;

% Określenie warunków początkowych:
y(1:kk)=Ypp;
u(1:kk)=Upp;
yzad(1:kk)=-0.3;
yzad(150:kk)=10;
yzad(300:kk)=0.5;
yzad(450:kk)=1.5;

% Ustawienia parametrów:
D = 600; 
N = 200; 
Nu = 20; 
lambda = 5;

% Pozyskanie odpowiedzi skokowej:
u(20:kk) = 1;
y(1:kk)=Ypp;
for k=7:kk
    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
end
S=y(21:kk);
S(kk-20:kk)=S(end);

% Wyliczenie macierzy M:
M=zeros(N,Nu);
for i=1:N
    for j=1:Nu
        if (i>=j)
            if(i-j+1<=D)
                M(i,j)=S(i-j+1);
            else
                M(i,j)=S(D);
            end
        end
    end
end

% Wyliczenie macierzy Mp:
Mp = zeros(N, D-1);
for i=1:N
    for j=1:D-1
        if i + j <= D
            Mp(i, j) = S(i + j) - S(j);
        else
            Mp(i, j) = S(D) - S(j);
        end
    end
end

% Wyliczenie macierzy K:
YY=eye(N,N);
L = lambda * eye(Nu);
K = (M'*YY*M + L)^(-1) * M'*YY;

% Inicjalizacja wektora dUp:
dUp(1:D-1)=0;
dUp=dUp';
u(1:kk) = Upp;
y(1:kk) = Ypp;

% Algorytm DMC:
for k=7:kk
    % Pomiar wyjścia:
    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));
    % Wyliczenie wektora dUp:
    for i=1:D-1
        dUp(i)=0;
        if k-i-1>0
            dUp(i)=u(k-i)-u(k-1-i);
        end
    end
    % Wyznaczenie wektorów:
    Y = y(k) * ones(N,1);
    Yzad = yzad(k) * ones(N, 1);
    % Wyznaczenie przyrostu sterowania:
    dU = K * (Yzad - Y - Mp * dUp);
    % Wyznaczenie sterowania:
    u(k) = u(k-1) + dU(1);
    % Uzwględnienie ograniczeń:
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
%print('pust5dmc_u_ga.png','-dpng','-r400');

figure;
hold on
stairs((1:kk),y);
stairs((1:kk),yzad, ':');
title("yzad, y, E="+ num2str(E)); 
xlabel('k'); 
ylabel('y');
legend('y','yzad');
title("y, E="+ num2str(E))
%print('pust5dmc_y_ga.png','-dpng','-r400');