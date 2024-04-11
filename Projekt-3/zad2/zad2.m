%Ustawienia
clear;
kk=200;

Upp=0;
Ypp=0;

y(1:kk)=Ypp;
u(1:kk)=Upp;

U=[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];

Y=zeros(length(U),kk);
U_skok=zeros(length(U),kk);
Chr_stat_U = (-1:0.01:1);
Chr_stat(1:length(U))=0;

%Odpowiedzi skokowe
for i=1:length(U)
    u(20:kk)=U(i);
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));
    end
    Chr_stat(i)=y(end);
    Y(i,:)=y;
    U_skok(i,:)=u;
end

figure
hold on
for i=1:length(U)
    stairs(1:kk,Y(i,:))
end
title('Odpowiedzi skokowe procesu:')
xlabel('k')
ylabel('y')
legend('u=-1','u=-0,75','u=-0,5','u=-0,25', ...
    'u=0','u=0,25','u=0,5','u=0,75','u=1')
%print('zad2_y.png','-dpng','-r400')

figure
hold on
for i=1:length(U)
    stairs(1:kk,U_skok(i,:))
end
title('Skoki sterowania:')
xlabel('k')
ylabel('y')
%print('zad2_u.png','-dpng','-r400')

%Wyznaczenie charakterystyki statycznej
for i=1:length(Chr_stat_U)
    u(20:kk)=Chr_stat_U(i);
    for k=7:kk
        y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));
    end
    Chr_stat(i)=y(end);
end

figure
plot(Chr_stat_U,Chr_stat)
title('Charakterystyka statyczna procesu:')
xlabel('u')
ylabel('y')
%print('zad2_chr_stat.png','-dpng','-r400')