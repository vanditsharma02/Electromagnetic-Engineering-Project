%fprintf("2");
function [capacitance]= cap(N,A,B,gaussquad)

NA = N;
NB = N;

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

capacitance = 0;
for i = 1:plates
    capacitance = capacitance+Q(i);
end

