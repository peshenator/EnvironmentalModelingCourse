function Tnew = TemperatureConvectionDiffusion(u,v,T)

global imax jmax dx dy dt lambda TBC xc;
Tnew = T;

% SPACE LOOP
for i=1:imax
    for j=1:jmax

        % x flux
        if(i==1)
            % convective part of the flux
            fp = 0.5*u(i+1,j)*(T(i+1,j) + T(i  ,j)) - 0.5*abs(u(i+1,j))*(T(i+1,j)-T(i  ,j));
            fm = 0;
            % diffusive part of the flux
%             fp = fp - lambda*(T(i+1,j) - T(i  ,j))/dx;
        elseif(i==imax)
            % convective part of the flux
            fp = 0;
            fm = 0.5*u(i  ,j)*(T(i  ,j) + T(i-1,j)) - 0.5*abs(u(i  ,j))*(T(i  ,j)-T(i-1,j));
            % diffusive part of the flux
%             fm = fm - lambda*(T(i  ,j) - T(i-1,j))/dx;
        else
            % convective part of the flux
            fp = 0.5*u(i+1,j)*(T(i+1,j) + T(i  ,j)) - 0.5*abs(u(i+1,j))*(T(i+1,j)-T(i  ,j));
            fm = 0.5*u(i  ,j)*(T(i  ,j) + T(i-1,j)) - 0.5*abs(u(i  ,j))*(T(i  ,j)-T(i-1,j));
            % diffusive part of the flux
%             fp = fp - lambda*(T(i+1,j) - T(i  ,j))/dx;
%             fm = fm - lambda*(T(i  ,j) - T(i-1,j))/dx;
        end

        % y flux
        % TBC at the bottom withexponent function
        if(j==1)
            Tbc = TBC*exp(-0.5*(xc(i)-0.5)^2/0.1^2);
            gp = 0.5*v(i,j+1)*(T(i,j+1) + T(i,j  )) - 0.5*abs(v(i,j+1))*(T(i,j+1)-T(i,j  ));
            gm =0;
%             gm =-lambda*(T(i,j) - Tbc)/(dy/2);
%             gp = gp - lambda*(T(i,j+1)-T(i,j  ))/dy;
        elseif(j==jmax)
            gp = 0;
            gm = 0.5*v(i,j  )*(T(i,j  ) + T(i,j-1)) - 0.5*abs(v(i,j  ))*(T(i,j  )-T(i,j-1));
%             gm = gm - lambda*(T(i,j  )-T(i,j-1))/dy;
        else
            gp = 0.5*v(i,j+1)*(T(i,j+1) + T(i,j  )) - 0.5*abs(v(i,j+1))*(T(i,j+1)-T(i,j  ));
            gm = 0.5*v(i,j  )*(T(i,j  ) + T(i,j-1)) - 0.5*abs(v(i,j  ))*(T(i,j  )-T(i,j-1));
%             gp = gp - lambda*(T(i,j+1)-T(i,j  ))/dy;
%             gm = gm - lambda*(T(i,j  )-T(i,j-1))/dy;
        end
        
        % finite volume update of the temperature
        Tnew(i,j) =  T(i,j) - dt/dx*( fp - fm ) - dt/dy*( gp - gm );
    end
end



end