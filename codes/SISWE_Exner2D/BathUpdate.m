
function bathbNew = BathUpdate(bathb,u,v,H)

    global dx dy dt Nx Ny;

    [fx,fy] = bflux(u,v,H);

    bathbNew = bathb - dt/dx*(fx(2:Nx+1,:)-fx(1:Nx,:)) - dt/dy*(fy(:,2:Ny+1)-fy(:,1:Ny));

end