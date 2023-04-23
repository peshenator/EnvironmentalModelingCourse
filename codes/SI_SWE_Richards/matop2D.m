function Apsi=matop2D(psi)
global dx dz dt K Kx Kz Nx Nz di H KL g
% diagonal part including the derivatives of the nonlinear functions 
Apsi = di.*psi;

fzm = zeros(Nz+1,Nx);
fzp = zeros(Nz+1,Nx);
fxm = zeros(Nz+1,Nx);
fxp = zeros(Nz+1,Nx);

% z fluxes
fzm(1     ,:) = 0;
fzm(2:Nz+1,:) = Kz(2:Nz+1 ,:).*(psi(2:Nz+1,:) - psi(1:Nz,:))/dz;
fzp(1:Nz  ,:) = fzm(2:Nz+1,:);
fzp(Nz+1,:) = 0;
% <-- Explanation: for the free surface layer, there is no mass flux across the
% free surface, but a mass flux across the bottom, given by
% Richards equation ! 

%  x flux
fxm(:,1     ) = 0;
fxm(:,2:Nx  ) = Kx(:,2:Nx).*(psi(:,2:Nx) - psi(:,1:Nx-1))/dx;
fxp(:,1:Nx-1) = fxm(:,2:Nx);
fxp(:,  Nx  ) = 0;

Apsi = Apsi - dt/dz*( fzp - fzm ) - dt/dx*( fxp - fxm );


% Apsi(1,   :) = Apsi(1   ,:) - dt/dz^2*( Kz(2     ,:).*(psi(2     ,:) - psi(1     ,:)) - 0 );
% Apsi(2:Nz,:) = Apsi(2:Nz,:) - dt/dz^2*( Kz(3:Nz+1,:).*(psi(3:Nz+1,:) - psi(2:Nz  ,:))...
%                                        -Kz(2:Nz  ,:).*(psi(2:Nz  ,:) - psi(1:Nz-1,:)) );  
% Apsi(Nz+1,:) = Apsi(Nz+1,:) - dt/dz^2*( 0 - Kz(Nz+1,:).*(psi(Nz+1,:) - psi(Nz,:)) );  

% Apsi(:,1) = Apsi(:,1) - dt/dx^2*( Kx(:,2).*(psi(:,2) - psi(:,1)) - 0 ); 
% Apsi(:,2:Nx-1) = Apsi(:,2:Nx-1) - dt/dx^2*( Kx(:,3:Nx).*(psi(:,3:Nx) - psi(:,2:Nx-1))...
%                                     -Kx(:,2:Nx-1).*(psi(:,2:Nx-1)-psi(:,1:Nx-2)) ); 
% Apsi(:,Nx) = Apsi(:,Nx) - dt/dx^2*( 0 - Kx(:,Nx).*(psi(:,Nx) - psi(:,Nx-1)) ); 


% % linear part 
% for i=1:Nx
%     for k=1:Nz+1
%         % z fluxes 
%         if(k==1)
%             Apsi(k,i) = Apsi(k,i)-dt/dz^2*( Kz(k+1,i)*(psi(k+1,i)-psi(k,i)) - 0 );  
%         elseif(k==Nz+1)
%             % for the free surface layer, there is no mass flux across the
%             % free surface, but a mass flux across the bottom, given by
%             % Richards equation ! 
%             Apsi(k,i) = Apsi(k,i)-dt/dz^2*( 0 - Kz(k,i)*(psi(k,i)-psi(k-1,i)) );  
%         else
%             Apsi(k,i) = Apsi(k,i)-dt/dz^2*( Kz(k+1,i)*(psi(k+1,i)-psi(k,i))...
%                                            -Kz(k,i)*(psi(k,i)-psi(k-1,i)) );  
%         end
%         % x fluxes for the Richards and the free surface equation 
%         if(i==1)
%             Apsi(k,i) = Apsi(k,i) - dt/dx^2*( Kx(k,i+1)*(psi(k,i+1)-psi(k,i)) - 0 ); 
%         elseif(i==Nx)
%             Apsi(k,i) = Apsi(k,i) - dt/dx^2*( 0 - Kx(k,i)*(psi(k,i)-psi(k,i-1)) ); 
%         else
%             Apsi(k,i) = Apsi(k,i) - dt/dx^2*( Kx(k,i+1)*(psi(k,i+1)-psi(k,i))...
%                                              -Kx(k,i)*(psi(k,i)-psi(k,i-1)) ); 
%         end
%     end
% end

end