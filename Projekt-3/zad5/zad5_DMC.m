% Inicjalizacja zmiennych:
clear;
Tp=0.5;

kk=600;
Upp=0;
Ypp=0;

Umax = 1;
Umin = -1;
E = 0;

% Inicjalizacja funkcji przynależności i pozyskanie
% odpowiedzi skokowych dla danej liczby regulatorów
% w wyznaczonych punktach pracy:
liczba_regulatorow=5;
szerokosc=Umax-Umin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if liczba_regulatorow==2
    N(1:liczba_regulatorow)=[35 10];
    Nu(1:liczba_regulatorow) = [1 1];              
    lambda(1:liczba_regulatorow) = [1 1];
    bok = szerokosc/(liczba_regulatorow+1);
    mf = [fismf("trapmf",[-100 -1 -1+bok -1+2*bok]);
      fismf("trapmf",[-1+bok -1+bok*2 1 100]);
      ];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=-0.66;
    y(1:kk)=-0.283777;
    u(20:kk) =-0.66+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S1=(y(21:kk)-y(1))/(u(end)-u(1));
    S1(kk-20:kk)=S1(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=0.66;
    y(1:kk)=5.82982;
    u(20:kk) =0.66+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S2=(y(21:kk)-y(1))/(u(end)-u(1));
    S2(kk-20:kk)=S2(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if liczba_regulatorow==3
    N(1:liczba_regulatorow)=[45 10 10];
    Nu(1:liczba_regulatorow) = [1 1 1];              
    lambda(1:liczba_regulatorow) = [1 1 1];
    bok = szerokosc/(liczba_regulatorow+1);
    mf = [fismf("trapmf",[-100 -1 -1+bok -1+2*bok]);
      fismf("trapmf",[-1+bok -1+bok*2 -1+bok*2 -1+bok*3]);
      fismf("trapmf",[-1+bok*2 -1+bok*3 1 100]);];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=-0.75;
    y(1:kk)=-0.291798;
    u(20:kk) =-0.75+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S1=(y(21:kk)-y(1))/(u(end)-u(1));
    S1(kk-20:kk)=S1(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=0;
    y(1:kk)=0;
    u(20:kk) =0+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S2=(y(21:kk)-y(1))/(u(end)-u(1));
    S2(kk-20:kk)=S2(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=0.75;
    y(1:kk)=7.28582;
    u(20:kk) =0.75+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S3=(y(21:kk)-y(1))/(u(end)-u(1));
    S3(kk-20:kk)=S3(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if liczba_regulatorow==4
    N(1:liczba_regulatorow)=[35 6 10 10];
    Nu(1:liczba_regulatorow) = [1 1 1 1];              
    lambda(1:liczba_regulatorow) = [1 1 1 1];
    bok = szerokosc/(liczba_regulatorow+1);
    mf = [fismf("trapmf",[-100 -1 -1+bok -1+2*bok]);
      fismf("trapmf",[-1+bok -1+bok*2 -1+bok*2 -1+bok*3]);
      fismf("trapmf",[-1+bok*2 -1+bok*3 -1+bok*3 -1+bok*4]);
      fismf("trapmf",[-1+bok*3 -1+bok*4 1 100]);
];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=-0.8;
    y(1:kk)=-0.296365;
    u(20:kk) =-0.8+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S1=(y(21:kk)-y(1))/(u(end)-u(1));
    S1(kk-20:kk)=S1(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=-0.2;
    y(1:kk)=-0.212438;
    u(20:kk) =-0.2+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S2=(y(21:kk)-y(1))/(u(end)-u(1));
    S2(kk-20:kk)=S2(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=0.2;
    y(1:kk)=0.733765;
    u(20:kk) =0.2+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S3=(y(21:kk)-y(1))/(u(end)-u(1));
    S3(kk-20:kk)=S3(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=0.8;
    y(1:kk)=8.14098;
    u(20:kk) =0.8+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S4=(y(21:kk)-y(1))/(u(end)-u(1));
    S4(kk-20:kk)=S4(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if liczba_regulatorow==5
    N(1:liczba_regulatorow)=[80 10 10 10 10];
    Nu(1:liczba_regulatorow) = [1 1 1 1 1];              
    lambda(1:liczba_regulatorow) = [1 1 1 1 1];
    bok = szerokosc/(liczba_regulatorow+1);
    mf = [fismf("trapmf",[-100 -1 -1+bok -1+2*bok]);
      fismf("trapmf",[-1+bok -1+bok*2 -1+bok*2 -1+bok*3]);
      fismf("trapmf",[-1+bok*2 -1+bok*3 -1+bok*3 -1+bok*4]);
      fismf("trapmf",[-1+bok*3 -1+bok*4 -1+bok*4 -1+bok*5]);
      fismf("trapmf",[-1+bok*4 -1+bok*5 1 100]);
      ];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=-0.833;
    y(1:kk)=-0.299143;
    u(20:kk) =-0.834+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S1=(y(21:kk)-y(1))/(u(end)-u(1));
    S1(kk-20:kk)=S1(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=-0.33;
    y(1:kk)=-0.249782;
    u(20:kk) =-0.33+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S2=(y(21:kk)-y(1))/(u(end)-u(1));
    S2(kk-20:kk)=S2(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=0;
    y(1:kk)=0;
    u(20:kk) =0+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S3=(y(21:kk)-y(1))/(u(end)-u(1));
    S3(kk-20:kk)=S3(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=0.33;
    y(1:kk)=1.68813;
    u(20:kk) =0.333+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S4=(y(21:kk)-y(1))/(u(end)-u(1));
    S4(kk-20:kk)=S4(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u(1:kk)=0.83;
    y(1:kk)=8.66867;
    u(20:kk) =0.834+0.01;
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
    end
    S5=(y(21:kk)-y(1))/(u(end)-u(1));
    S5(kk-20:kk)=S5(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
x = -1:0.001:1;
mi = evalmf(mf,x);
D=600;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wyliczenie macierzy M:
if liczba_regulatorow==2
    M1=zeros(N(1),Nu(1));
    M2=zeros(N(2),Nu(2));
    for i=1:N(1)
        for j=1:Nu(1)
            if (i>=j)
                if(i-j+1<=D)
                    M1(i,j)=S1(i-j+1);
                else
                    M1(i,j)=S1(D);
                end
            end
        end
    end
    for i=1:N(2)
        for j=1:Nu(2)
            if (i>=j)
                if(i-j+1<=D)
                    M2(i,j)=S2(i-j+1);
                else
                    M2(i,j)=S2(D);
                end
            end
        end
    end
end
if liczba_regulatorow==3
    M1=zeros(N(1),Nu(1));
    M2=zeros(N(2),Nu(2));
    M3=zeros(N(3),Nu(3));
    for i=1:N(1)
        for j=1:Nu(1)
            if (i>=j)
                if(i-j+1<=D)
                    M1(i,j)=S1(i-j+1);
                else
                    M1(i,j)=S1(D);
                end
            end
        end
    end
    for i=1:N(2)
        for j=1:Nu(2)
            if (i>=j)
                if(i-j+1<=D)
                    M2(i,j)=S2(i-j+1);
                else
                    M2(i,j)=S2(D);
                end
            end
        end
    end
    for i=1:N(3)
        for j=1:Nu(3)
            if (i>=j)
                if(i-j+1<=D)
                    M3(i,j)=S3(i-j+1);
                else
                    M3(i,j)=S3(D);
                end
            end
        end
    end
end
if liczba_regulatorow==4
    M1=zeros(N(1),Nu(1));
    M2=zeros(N(2),Nu(2));
    M3=zeros(N(3),Nu(3));
    M4=zeros(N(4),Nu(4));
    for i=1:N(1)
        for j=1:Nu(1)
            if (i>=j)
                if(i-j+1<=D)
                    M1(i,j)=S1(i-j+1);
                else
                    M1(i,j)=S1(D);
                end
            end
        end
    end
    for i=1:N(2)
        for j=1:Nu(2)
            if (i>=j)
                if(i-j+1<=D)
                    M2(i,j)=S2(i-j+1);
                else
                    M2(i,j)=S2(D);
                end
            end
        end
    end
    for i=1:N(3)
        for j=1:Nu(3)
            if (i>=j)
                if(i-j+1<=D)
                    M3(i,j)=S3(i-j+1);
                else
                    M3(i,j)=S3(D);
                end
            end
        end
    end
    for i=1:N(4)
        for j=1:Nu(4)
            if (i>=j)
                if(i-j+1<=D)
                    M4(i,j)=S4(i-j+1);
                else
                    M4(i,j)=S4(D);
                end
            end
        end
    end
end
if liczba_regulatorow==5
    M1=zeros(N(1),Nu(1));
    M2=zeros(N(2),Nu(2));
    M3=zeros(N(3),Nu(3));
    M4=zeros(N(4),Nu(4));
    M5=zeros(N(5),Nu(5));
    for i=1:N(1)
        for j=1:Nu(1)
            if (i>=j)
                if(i-j+1<=D)
                    M1(i,j)=S1(i-j+1);
                else
                    M1(i,j)=S1(D);
                end
            end
        end
    end
    for i=1:N(2)
        for j=1:Nu(2)
            if (i>=j)
                if(i-j+1<=D)
                    M2(i,j)=S2(i-j+1);
                else
                    M2(i,j)=S2(D);
                end
            end
        end
    end
    for i=1:N(3)
        for j=1:Nu(3)
            if (i>=j)
                if(i-j+1<=D)
                    M3(i,j)=S3(i-j+1);
                else
                    M3(i,j)=S3(D);
                end
            end
        end
    end
    for i=1:N(4)
        for j=1:Nu(4)
            if (i>=j)
                if(i-j+1<=D)
                    M4(i,j)=S4(i-j+1);
                else
                    M4(i,j)=S4(D);
                end
            end
        end
    end
    for i=1:N(5)
        for j=1:Nu(5)
            if (i>=j)
                if(i-j+1<=D)
                    M5(i,j)=S5(i-j+1);
                else
                    M5(i,j)=S5(D);
                end
            end
        end
    end
end

% Wyliczenie macierzy Mp:
if liczba_regulatorow==2
    Mp1 = zeros(N(1), D-1);
    Mp2 = zeros(N(2), D-1);
    for i=1:N(1)
        for j=1:D-1
            if i + j <= D
                Mp1(i, j) = S1(i + j) - S1(j);
            else
                Mp1(i, j) = S1(D) - S1(j);
            end
        end
    end
    for i=1:N(2)
        for j=1:D-1
            if i + j <= D
                Mp2(i, j) = S2(i + j) - S2(j);
            else
                Mp2(i, j) = S2(D) - S2(j);
            end
        end
    end
end
if liczba_regulatorow==3
    Mp1 = zeros(N(1), D-1);
    Mp2 = zeros(N(2), D-1);
    Mp3 = zeros(N(3), D-1);
    for i=1:N(1)
        for j=1:D-1
            if i + j <= D
                Mp1(i, j) = S1(i + j) - S1(j);
            else
                Mp1(i, j) = S1(D) - S1(j);
            end
        end
    end
    for i=1:N(2)
        for j=1:D-1
            if i + j <= D
                Mp2(i, j) = S2(i + j) - S2(j);
            else
                Mp2(i, j) = S2(D) - S2(j);
            end
        end
    end
    for i=1:N(3)
        for j=1:D-1
            if i + j <= D
                Mp3(i, j) = S3(i + j) - S3(j);
            else
                Mp3(i, j) = S3(D) - S3(j);
            end
        end
    end
end
if liczba_regulatorow==4
    Mp1 = zeros(N(1), D-1);
    Mp2 = zeros(N(2), D-1);
    Mp3 = zeros(N(3), D-1);
    Mp4 = zeros(N(4), D-1);
    for i=1:N(1)
        for j=1:D-1
            if i + j <= D
                Mp1(i, j) = S1(i + j) - S1(j);
            else
                Mp1(i, j) = S1(D) - S1(j);
            end
        end
    end
    for i=1:N(2)
        for j=1:D-1
            if i + j <= D
                Mp2(i, j) = S2(i + j) - S2(j);
            else
                Mp2(i, j) = S2(D) - S2(j);
            end
        end
    end
    for i=1:N(3)
        for j=1:D-1
            if i + j <= D
                Mp3(i, j) = S3(i + j) - S3(j);
            else
                Mp3(i, j) = S3(D) - S3(j);
            end
        end
    end
    for i=1:N(4)
        for j=1:D-1
            if i + j <= D
                Mp4(i, j) = S4(i + j) - S4(j);
            else
                Mp4(i, j) = S4(D) - S4(j);
            end
        end
    end
end
if liczba_regulatorow==5
    Mp1 = zeros(N(1), D-1);
    Mp2 = zeros(N(2), D-1);
    Mp3 = zeros(N(3), D-1);
    Mp4 = zeros(N(4), D-1);
    Mp5 = zeros(N(5), D-1);
    for i=1:N(1)
        for j=1:D-1
            if i + j <= D
                Mp1(i, j) = S1(i + j) - S1(j);
            else
                Mp1(i, j) = S1(D) - S1(j);
            end
        end
    end
    for i=1:N(2)
        for j=1:D-1
            if i + j <= D
                Mp2(i, j) = S2(i + j) - S2(j);
            else
                Mp2(i, j) = S2(D) - S2(j);
            end
        end
    end
    for i=1:N(3)
        for j=1:D-1
            if i + j <= D
                Mp3(i, j) = S3(i + j) - S3(j);
            else
                Mp3(i, j) = S3(D) - S3(j);
            end
        end
    end
    for i=1:N(4)
        for j=1:D-1
            if i + j <= D
                Mp4(i, j) = S4(i + j) - S4(j);
            else
                Mp4(i, j) = S4(D) - S4(j);
            end
        end
    end
    for i=1:N(5)
        for j=1:D-1
            if i + j <= D
                Mp5(i, j) = S5(i + j) - S5(j);
            else
                Mp5(i, j) = S5(D) - S5(j);
            end
        end
    end
end

% Wyliczenie macierzy K:
if liczba_regulatorow==2
    K1 = (M1'*M1 + lambda(1) * eye(Nu(1)))^(-1) * M1';
    K2 = (M2'*M2 + lambda(2) * eye(Nu(2)))^(-1) * M2';
end
if liczba_regulatorow==3
    K1 = (M1'*M1 + lambda(1) * eye(Nu(1)))^(-1) * M1';
    K2 = (M2'*M2 + lambda(2) * eye(Nu(2)))^(-1) * M2';
    K3 = (M3'*M3 + lambda(3) * eye(Nu(3)))^(-1) * M3';
end
if liczba_regulatorow==4
    K1 = (M1'*M1 + lambda(1) * eye(Nu(1)))^(-1) * M1';
    K2 = (M2'*M2 + lambda(2) * eye(Nu(2)))^(-1) * M2';
    K3 = (M3'*M3 + lambda(3) * eye(Nu(3)))^(-1) * M3';
    K4 = (M4'*M4 + lambda(4) * eye(Nu(4)))^(-1) * M4';
end
if liczba_regulatorow==5
    K1 = (M1'*M1 + lambda(1) * eye(Nu(1)))^(-1) * M1';
    K2 = (M2'*M2 + lambda(2) * eye(Nu(2)))^(-1) * M2';
    K3 = (M3'*M3 + lambda(3) * eye(Nu(3)))^(-1) * M3';
    K4 = (M4'*M4 + lambda(4) * eye(Nu(4)))^(-1) * M4';
    K5 = (M5'*M5 + lambda(5) * eye(Nu(5)))^(-1) * M5';
end

% Inicjalizacja wektora dUp:
dUp(1:D-1)=0;
dUp=dUp';
% Określenie warunków początkowych:
y(1:kk)=Ypp;
u(1:kk)=Upp;
yzad(1:kk)=-0.3;
yzad(150:kk)=10;
yzad(300:kk)=0.5;
yzad(450:kk)=1.5;

% Algorytm DMC:
for k=7:kk
    % Pomiar wyjścia:
    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));
    % Wyliczenie stopnia przynależności:
    mi=evalmf(mf,u(k-1));
    % Wyliczenie wektora dUp:
    for i=1:D-1
        dUp(i)=0;
        if k-i-1>0
            dUp(i)=u(k-i)-u(k-1-i);
        end
    end
    % Wyliczenie przyrostu sterowania:
    if liczba_regulatorow==2
        du1=K1*(yzad(k)*ones(N(1),1)-y(k)*ones(N(1),1)-Mp1*dUp);
        du2=K2*(yzad(k)*ones(N(2),1)-y(k)*ones(N(2),1)-Mp2*dUp);
        dU = (du1(1)*mi(1)+du2(1)*mi(2))/(mi(1)+mi(2));
    end
    if liczba_regulatorow==3
        du1=K1*(yzad(k)*ones(N(1),1)-y(k)*ones(N(1),1)-Mp1*dUp);
        du2=K2*(yzad(k)*ones(N(2),1)-y(k)*ones(N(2),1)-Mp2*dUp);
        du3=K3*(yzad(k)*ones(N(3),1)-y(k)*ones(N(3),1)-Mp3*dUp);
        dU = (du1(1)*mi(1)+du2(1)*mi(2)+ ...
            du3(1)*mi(3))/(mi(1)+mi(2)+mi(3));
    end
    if liczba_regulatorow==4
        du1=K1*(yzad(k)*ones(N(1),1)-y(k)*ones(N(1),1)-Mp1*dUp);
        du2=K2*(yzad(k)*ones(N(2),1)-y(k)*ones(N(2),1)-Mp2*dUp);
        du3=K3*(yzad(k)*ones(N(3),1)-y(k)*ones(N(3),1)-Mp3*dUp);
        du4=K4*(yzad(k)*ones(N(4),1)-y(k)*ones(N(4),1)-Mp4*dUp);
        dU = (du1(1)*mi(1)+du2(1)*mi(2)+du3(1)*mi(3)+ ...
            du4(1)*mi(4))/(mi(1)+mi(2)+mi(3)+mi(4));
    end
    if liczba_regulatorow==5
        du1=K1*(yzad(k)*ones(N(1),1)-y(k)*ones(N(1),1)-Mp1*dUp);
        du2=K2*(yzad(k)*ones(N(2),1)-y(k)*ones(N(2),1)-Mp2*dUp);
        du3=K3*(yzad(k)*ones(N(3),1)-y(k)*ones(N(3),1)-Mp3*dUp);
        du4=K4*(yzad(k)*ones(N(4),1)-y(k)*ones(N(4),1)-Mp4*dUp);
        du5=K5*(yzad(k)*ones(N(5),1)-y(k)*ones(N(5),1)-Mp5*dUp);
        dU = (du1(1)*mi(1)+du2(1)*mi(2)+du3(1)*mi(3)+ ...
            du4(1)*mi(4)+du5(1)*mi(5))/(mi(1)+ ...
            mi(2)+mi(3)+mi(4)+mi(5));
    end
    % Wyliczenie sterowania:
    u(k) = u(k-1) + dU;
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
print('pust_dmc_u_5_regulatory_przed.png','-dpng','-r400');

figure;
hold on
stairs((1:kk),y);
stairs((1:kk),yzad, ':');
title("yzad, y, E="+ num2str(E)); 
xlabel('k'); 
ylabel('y');
legend('y','yzad');
title("y, E="+ num2str(E))
print('pust_dmc_y_5_regulatory_przed.png','-dpng','-r400');