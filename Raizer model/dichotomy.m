function[answer, delta, iternum ] = dichotomy(Z,T,V,x, L,R,eps)

M = 100;
a = zeros(1,100);
b = zeros(1,100);
c = zeros(1,100);

m = 1;
a(1) = L;
b(1) = R;

if equation(Z,T,V,x, a(1)) == 0 % ��������, �� �������� �� ��������� ����� �������
    answer = a(1);
    delta = 0;
elseif equation(Z,T,V,x, b(1)) == 0
    answer = b(1);
    delta = 0;
else

    delta = b(1)-a(1);

    ab = equation( Z,T,V,x, a(m) )*equation( Z,T,V,x, b(m) );
    if ab > 0 % ���� � ��������� ������ �������� ������ �����, �� ������ ���� ���, ���� ���������
        disp('Error. No roots or multiple roots.');
        answer = -1;
    else % ���� � ���� ����, �� �������� ��������

        while delta > eps
            c(m) = 0.5*a(m) + 0.5*b(m); % �������� �������
            if equation(Z,T,V,x, a(m)) == 0 % �� ������ �������� ���������, �� �������� �� ������� � �������� �������
                answer = a(m);
                delta = 0;
                break;
            elseif equation(Z,T,V,x, b(m)) == 0
                answer = b(m);
                delta = 0;
                break;
            elseif equation(Z,T,V,x, c(m)) == 0
                answer = c(m);
                delta = 0;
                break;
            end

            ac = equation( Z,T,V,x, a(m) )*equation( Z,T,V,x, c(m) );
            cb = equation( Z,T,V,x, b(m) )*equation( Z,T,V,x, c(m) );
            if ac < 0
                b(m+1) = c(m);
                a(m+1) = a(m);
            elseif cb < 0
                a(m+1) = c(m);
                b(m+1) = b(m);
            end
            answer = 0.5*a(m+1) + 0.5*b(m+1); % ��� ��������� � �������� ������, ����� ������ �� �����
            delta = b(m+1) - a(m+1);
            m = m+1;
            if( m >= M ) % �� ������������
                break;
            end
        end
    end
end
iternum = m-1;
end