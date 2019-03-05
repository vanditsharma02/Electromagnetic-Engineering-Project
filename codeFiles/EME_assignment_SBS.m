A = 1;
B = 1;

NA = 50;
NB = 50;

gaussquad = 5;
H = 2;

eps = 8.85*(10^(-12));
pie = 3.14;
plates = NA*NB;

P = zeros(plates*2,plates*2);

fun = @(x,y) 1./(sqrt(x.^2 + y.^2));
q = integral2(fun,-A/(NA*2),A/(NA*2),-B/(NB*2),B/(NB*2));




%disp(q);
%temp = 8*log(1+(2^0.5))*(A/(2*NA));

X1 = [];
Y1 = [];
X2 = [];
Y2 = [];

[u,v] = lgwt(gaussquad,-(A/(2*NA)),(A/(2*NA)));
disp(u);
disp(v);

for i = 1:(plates)
for j = 1:(plates*2)
		
        if j == i
            P(i,j) =(1/((4*pie*eps)))*q;

        else
            if(j<=plates)
                temp = 0;
                for k = 1:gaussquad
                    for l = 1:gaussquad
                        x1 = [((mod(j-1,NA)+1/2)*(B/NB)+u(k)) ((floor((j-1)/NA)+1/2)*(A/NA)+u(l))];
                        x2 = [((mod(i-1,NA)+1/2)*(B/NB)) ((floor((i-1)/NA)+1/2)*(A/NA))];
                        dist = (((x1(1)-x2(1))^2)+((x1(2)-x2(2))^2))^0.5;
                        temp = temp+ v(k)*v(l)*(1/((4*pie*eps)*(dist)));
                    end
                end
                P(i,j) = temp;
            else
                temp = 0;
                for k = 1:gaussquad
                    for l = 1:gaussquad
                        x1 = [((mod(j-1-plates,NA)+1/2)*(B/NB)+u(k)+H) ((floor((j-1-plates)/NA)+1/2)*(A/NA)+u(l))];
                        x2 = [((mod(i-1,NA)+1/2)*(B/NB)) ((floor((i-1)/NA)+1/2)*(A/NA))];
                        dist = (((x1(1)-x2(1))^2)+((x1(2)-x2(2))^2)+(H^2))^0.5;
                        temp = temp+ v(k)*v(l)*(1/((4*pie*eps)*(dist)));
                    end
                end
                P(i,j) = temp;
            end    
        end
    end
    X1 = [X1,(mod(i-1,NA)+1/2)*(B/NB)];
    Y1 = [Y1,(floor((i-1)/NA)+1/2)*(A/NA)];
end

for i= plates+1:plates*2
    for j = plates+1: plates*2
        P(i,j) = P(i-plates,j-plates);
    end
end

for i= plates+1:plates*2
    for j = 1: plates
        P(i,j) = P(i-plates,j);
    end
end

for i=1:NA*NB
    X2(i) = X1(i)+H;
end
Y2 = Y1;
%disp(P);
lin=[];
for j = 1:plates
    lin(j,1)=1/2;
    lin(j+plates,1)=-1/2;
end
%disp(lin);
P = P./((B/NB)*(A/NA));
Q = inv(P)*lin;
Q1 = [];
Q2 = [];
for j = 1:plates
    Q1(j) = Q(j);
    Q2(j) = Q(j+plates);
end

surfQ1=[,];
for i = 1:NA
	for j = 1:NB
        surfQ1(i,j) = Q1(((i-1)*NA)+j);
    end
end       

surfQ2=[,];
for i = 1:NA
	for j = 1:NB
        surfQ2(i,j) = Q2(((i-1)*NA)+j);
    end
end     
        
disp(Q1);
plot3(X1,Y1,Q1);

hold on;
plot3(X2,Y2,Q2);
hold off;

figure;
surf(surfQ1);

hold on;
surf(surfQ2);
hold off;
%%%%%%%%%%%%%%%%%%%%%% PART 2-FINDING ELECTRIC FIELD %%%%%%%%%%%%%%%%%%%%%%

nx = 100;
ny = 100;

for i = 1:nx
    for j = 1:ny
        s(i,j) = -1+ ((5*(i))/nx);
    end
end
    
for j = 1:ny
    for i = 1:nx
        t(i,j) = -1+ ((3*(j))/ny);
    end
end

V = zeros(nx,ny);
for j = 1:nx
    for k = 1:ny
        for i = 1:2*plates
            if i<= plates
                V(j,k) = V(j,k)+(1/(4*pie*eps))*Q(i)*((s(j,k)-X1(i))^2+(t(j,k)-Y1(i))^2+(0.001)^2)^(-0.5);
            else
                V(j,k) = V(j,k)+(1/(4*pie*eps))*Q(i)*((s(j,k)-X2(i-plates))^2+(t(j,k)-Y2(i-plates))^2+(0.001)^2)^(-0.5);
            end
        end
    end
end

figure;
plot3(s,t,V);


[c,d,e] = meshgrid(-1:0.2:2+H,-1:0.2:2,-2.2:0.4:2.2);
u=0;
for i = 1:plates
    u = u +(1/(4*pie*eps)).*(c-X1(i)).*Q1(i).*((c-X1(i)).^2+(d-Y1(i)).^2+(e).^2).^(-1.5);
    u = u +(1/(4*pie*eps)).*(c-X2(i)).*Q2(i).*((c-X2(i)).^2+(d-Y2(i)).^2+(e).^2).^(-1.5);
end
v=0;
for i = 1:plates
    v = v +(1/(4*pie*eps)).*(d-Y1(i)).*Q1(i).*((c-X1(i)).^2+(d-Y1(i)).^2+(e).^2).^(-1.5);
    v = v +(1/(4*pie*eps)).*(d-Y2(i)).*Q2(i).*((c-X2(i)).^2+(d-Y2(i)).^2+(e).^2).^(-1.5);
end
w=0;
for i = 1:plates
    w = w +(1/(4*pie*eps)).*(e).*Q1(i).*((c-X1(i)).^2+(d-Y1(i)).^2+(e).^2).^(-1.5);
    w = w +(1/(4*pie*eps)).*(e).*Q2(i).*((c-X2(i)).^2+(d-Y2(i)).^2+(e).^2).^(-1.5); 
end

figure;
quiver3(c,d,e,u,v,w,'linewidth',1);

%figure;
%stream3(c,d,e,u,v,w,startx,starty,startz);

capacitance = 0;
for i = 1:plates
    capacitance = capacitance+Q1(i);
end
capacitance

