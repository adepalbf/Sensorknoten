Nx = 2;
Ny = 1;
dx = 1;
dy = 1;

N = Nx * Ny;

ULA = zeros(3,N);
for ii = 1:Nx
    for jj = 1:Ny
        ULA(1,jj+Ny*(ii-1)) = (ii-1)*dx;
        ULA(2,jj+Ny*(ii-1)) = (jj-1)*dy;
    end
end
            