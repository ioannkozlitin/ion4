function data(Z)
global PHI;
global M_EL;

db = 2; % 1 -> TEFIS, else -> NIST

el_num = length(Z);
Z_max = max(Z);

file_title = 'masses.txt';
fid = fopen(file_title,'r');
S1 = fscanf(fid,'%g');
fclose(fid);

if( db == 1 )
    file_title = 'TEFIS.txt';
    fid = fopen(file_title,'r');
    S2 = fscanf(fid,'%g');
    fclose(fid);

    phi0 = zeros(103);
    num = 0;
    for m = 1:103
        for n = 1:m
            phi0(m,n) = S2(num + n + 1);
        end
        num = num + m + 1;
    end
else
    file_title = 'NIST.txt';
    fid = fopen(file_title,'r');
    S2 = fscanf(fid,'%g');
    fclose(fid);

    phi0 = zeros(103);
    num = 0;
    for m = 1:103
        for n = 1:m
            phi0(m,n) = S2(num + n);
        end
        num = num + m;
    end
end


PHI  = zeros(Z_max,el_num);
M_EL = zeros(1,el_num);
for j = 1:el_num
    for n = 1:Z(j)
        PHI(n,j) = phi0( Z(j),n )/27.26;
    end
    M_EL(j) = S1( Z(j) );
end
end