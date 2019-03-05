%fprintf("2");

A = 1;
B = 1;

NA = 50;
NB = 50;

gaussquad = 3;


eps = 8.85*(10^(-12));
pie = 3.14;
plates = NA*NB;

P = zeros(plates,plates);

fun = @(x,y) 1./(sqrt(x.^2 + y.^2));
q = integral2(fun,-A/(NA*2),A/(NA*2),-B/(NB*2),B/(NB*2));
%disp(q);
%temp = 8*log(1+(2^0.5))*(A/(2*NA));

X = [];
Y = [];

[u,v] = lgwt(gaussquad,-(A/(2*NA)),(A/(2*NA)));
disp(u);
disp(v);

for i = 1:plates
	for j = 1:plates
		
        if j == i
            P(i,j) = (1/((4*pie*eps)))*q;

        else
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
        end
    end
    X = [X,(mod(i-1,NA)+1/2)*(B/NB)];
    Y = [Y,(floor((i-1)/NA)+1/2)*(A/NA)];
end

%disp(P);

lin = ones(plates,1);
%disp(lin);
P = P./((B/NB)*(A/NA));
Q = pinv(P)*lin;

surfQ=[,];
for i = 1:NA
	for j = 1:NB
        surfQ(i,j) = Q(((i-1)*NA)+j);
    end
end       
        
disp(Q);
plot3(X,Y,Q);

figure;
surf(surfQ);

%%%%%%%%%%%%%%%%%%%%%% PART 2-FINDING ELECTRIC FIELD %%%%%%%%%%%%%%%%%%%%%%

nx = 100;
ny = 100;

for i = 1:nx
    for j = 1:ny
        s(i,j) = -1+ ((3*(i))/nx);
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
        for i = 1:plates
            V(j,k) = V(j,k)+(1/(4*pie*eps))*Q(i)*((s(j,k)-X(i))^2+(t(j,k)-Y(i))^2+(0.001)^2)^(-0.5);
        end
    end
end

figure;
plot3(s,t,V);

[c,d,e] = meshgrid(-2:0.2:2,-2:0.2:2,-2.2:0.4:2.2);
u=0;
for i = 1:plates
    u = u +(1/(4*pie*eps)).*(c-X(i)).*Q(i).*((c-X(i)).^2+(d-Y(i)).^2+(e).^2).^(-1.5);
end
v=0;
for i = 1:plates
    v = v +(1/(4*pie*eps)).*(d-Y(i)).*Q(i).*((c-X(i)).^2+(d-Y(i)).^2+(e).^2).^(-1.5);
end
w=0;
for i = 1:plates
    w = w +(1/(4*pie*eps)).*(e).*Q(i).*((c-X(i)).^2+(d-Y(i)).^2+(e).^2).^(-1.5);
end


figure;
quiver3(c,d,e,u,v,w,'linewidth',1);

%figure;
%stream3(c,d,e,u,v,w,startx,starty,startz);


capacitance = 0;
for i = 1:plates
    capacitance = capacitance+Q(i);
end

disp(capacitance);
