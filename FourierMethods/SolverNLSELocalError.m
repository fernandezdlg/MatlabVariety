function [inst] = SolverNLSELocalError(l,plotProp,rc,deltaG)%,Lz)
%SolverNLSELocalError is a solver to the Nonlinear Schroedinger Equation
%(NLSE)
%
% This version also consists on an implementation of the Fourier Split
% Step method with a Local Error algorithm which varies the size of the Z
% steps in order to achieve stable relative errors on propagation with a
% reduced number of steps. For more information please consult 
% O. V. Sinkin, R. Holzlohner, J. Zweck and C. R. Menyuk, "Optimization of
% the split-step Fourier method in modeling optical-fiber communications 
% systems," in Journal of Lightwave Technology, vol. 21, no. 1, pp. 61-68,
% Jan. 2003, doi: 10.1109/JLT.2003.808628.
%
% To modify the equation, one modifies the N(u) term whenever the Potential
% Operator handle function is called. The coefficients for the NLSE can be
% changed directly.
%
%  Inputs: l - Lamda, constant of propagation phase evolution
%          plotProp - 1 for viewing propagation, 0 to not view it
%          rc - Stability threshold (if the error gets greater than rc,
%               stability is lost)
%          deltaG - Goal local error.
%  Outputs: inst - the distance upon which stability is lost.
% Last Update: January 18, 2021
% Juan Antonio Fernandez de la Garza -- juanfernandezdlg@gmail.com

%% Computational framework initialization
% Equation to solve: iDz(u) +  A*Dx2(u) + B*N(u)u = 0
A = 0.5; % Typically 0.5
B = 1.0; % Typically 1.0

Lx = 25; % X-Length of simulation

Nx = 2048; % 2^(integer) number to perform efficient FFTs
disp(['Nx =',num2str(Nx)])

X = linspace(-Lx/2,Lx/2,Nx)'; % X-definition
dx = X(2)-X(1); % X-spacing

%% Propagation parameters definition
Lz = 10; % Z-Length of simulation
dz = dx^2/20; % Empirical parameter for the simulation (for small errors)

% Z blocks to reduce memory requirements (only Nbloques steps are stored)
Nbloques = 101; % Number of blocks
disp(['Number of blocks: ',num2str(Nbloques)])
disp(['Initial Z-spacing dz=',num2str(dz)])

Z = zeros(1,Nbloques);

inst = Lz; % Assume full stability

%% Soliton problem definition
u = zeros(Nx,Nbloques); % Array of Lorentzian solitons recorded instants

% Soliton initialization
u(:,1) = sech(X); % INSERT FUNCTION HERE

% Potential function
V = ones(size(u(:,1))); % INSERT FUNCTION HERE

%% Momentum discretization, Nx can be even or odd
K = linspace(0,2*pi/dx,Nx+1)'; % k-Space adjusted to match Matlab's FFT 
                               % convention
K = K(1:Nx);
dk = K(2);
K = fftshift(K);
K(1:floor(Nx/2)) = K(1:floor(Nx/2)) -(Nx)*dk; 
K = ifftshift(K);

% % Kinetic Operator
Ok = @(h) exp(-A*1i.*h.*K.^2);
% Potential Operator
Ov = @(h,N) exp(1i*B.*h.*N).*exp(1i*h.*V);

%% Section for tracking norm (power) (Comment if not of interest)
Ptot = zeros(Nbloques,1);
Ptot(1) = sum(abs(u(:,1)).^2)*dx;

%% To compute the error ||u(z)-u(0)e^(i lambda z)||/||u(0)||
uNorm = sqrt(sum(abs(u(:,1)).^2).*dx); % u norm
InstError = 0; % Instantaneous error variable definition

%% Section for introduction of gaussian noise 
% u(:,1) = awgn(u(:,1),50);
% shft = 1;
% u(shft:end,1) = u(1:end+1-shft,1);

%% Fourier split step propagation with Local-Error control via dz
% u([1,end],1) = [0;0]; % to set boundary conditions
ucoarse = u(:,1);
ufine = ucoarse;
u4 = ufine;
for bz = 2:Nbloques
    Z(bz) = Z(bz-1);
    while Z(bz) < (bz-1)*Lz/Nbloques
        accept = 0;
        while accept == 0
            ucoarse = u4;
            ufine = ucoarse;            
            %% Coarse propagation
            % Kinetic halfstep
            ucoarse = ifft(Ok(dz).*fft(ucoarse));
            % Potential step
            ucoarse = Ov(2*dz,abs(ucoarse).^2).*ucoarse; 
            % Kinetic halfstep
            ucoarse = ifft(Ok(dz).*fft(ucoarse));

            %% Fine propagation
            ufine = ifft(Ok(dz/2).*fft(ufine));
            % Potential step
            ufine = Ov(dz,abs(ufine).^2).*ufine; 
            % Kinetic fullstep
            ufine = ifft(Ok(dz).*fft(ufine));
            % Potential step
            ufine = Ov(dz,abs(ufine).^2).*ufine; 
            % Kinetic halfstep
            ufine = ifft(Ok(dz/2).*fft(ufine));
        
            %% Error check
            delta = sqrt(sum(abs(ufine-ucoarse).^2).*dx)/uNorm;
            
            if delta > 2*deltaG
                accept = 0;
                dz = dz/2;
            elseif delta > deltaG
                accept = 1;
                %% "Richardson term"
                u4 = (4*ufine-ucoarse)./3;
%                 u4([1,end]) = [0;0];
                Z(bz) = Z(bz) + 2*dz;
                dz = dz/(2^(1/3));
            elseif delta < deltaG/2
                accept = 1;
                %% "Richardson term"
                u4 = (4*ufine-ucoarse)./3;
%                 u4([1,end]) = [0;0];
                Z(bz) = Z(bz) + 2*dz;
                dz = dz*2^(1/3);
            else
                accept = 1;
                %% "Richardson term"
                u4 = (4*ufine-ucoarse)./3;
%                 u4([1,end]) = [0;0];
                Z(bz) = Z(bz) + 2*dz;
            end
        end
    end
    % Store propagation
    u(:,bz) = u4;
    
    % Power is tracked (Comment if not of interest)
    Ptot(bz) = sum(abs(u4).^2)*dx;
   
    % Check for instability:
    InstError = sqrt(sum(abs(u(:,1)*exp(1i*l*Z(bz))-u(:,bz)).^2).*dx)/uNorm;
    % disp(InstError)
    disp(dz)
    if InstError > rc
        inst = Z(bz);
        disp('Instability')
        break
    end 
end


%% Graphical results
if plotProp==1
    % Initial profile
    figure(1)
    plot(X,u(:,1))
    title('Initial soliton''s transverse amplitude')
    xlabel('$x$')
    ylabel('$u(x)$')
    axis([min(X) max(X) min(u(:,1))-range(u(:,1))*.1 max(u(:,1))+range(u(:,1))*.1])% Propagation

    % Potential function profile
    figure(2)
    plot(X,V)
    title('Potential function')
    xlabel('$x$')
    ylabel('$V(x)$')

    % Soliton propagation (abs value)
    figure(3)
    [XX,ZZ] = ndgrid(X,Z);
    s=surf(XX,ZZ,abs(u(:,:)));
    colormap(hot)
    s.EdgeColor = 'none'; 
    title('$|U(x,z)|$')
    view(2)
    xlim([min(X),max(X)])
    ylim([min(Z),max(Z)])
    xlabel('$x$')
    ylabel('$z$')
    shading interp
    
    % Power tracking
    figure(4)
    plot(Z,Ptot)
    title('Transverse power according to the propagation distance $P(z)$')
    axis([min(Z) max(Z) 0 max(Ptot)*1.1])
    xlabel('$z$')
    ylabel('$P(z)$')
end


end
