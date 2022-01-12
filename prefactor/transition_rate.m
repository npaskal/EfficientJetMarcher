function r_ab = transition_rate(flag,epsilon,a)
% if flag == 0, compute the transition rate using Transition Path Theory
% if flag == 1, compute transition rate using Bouchet-Reygner formula
%clear all
%close all
%epsilon = 0.01;
fsz = 24;
%flag = 0;
fig = 'y'; % fig = 'y' -- draw figures
%% Maier-Stein vector field
%a = 3;
MSx = @(x,y)x-x.^3-a*x.*y.^2;
MSy = @(x,y)-(1+x.^2).*y;
Pot = @(x,y)x.^4/4-x.^2/2+x.^2.*y.^2/2+y.^2/2;
MSdivb = @(x,y)-4*x.^2-a*y.^2;
if a == 10
    MAPdata = load('MSmap.mat');
    map = MAPdata.x;
end
%% compute the transition rate using the transition path theory
if flag == 0
    N = 1025;
    xmin = -1.5; xmax = 1.5;
    ymin = -0.75; ymax = 0.75;
    gx = linspace(xmin,xmax,N)';
    gy = linspace(ymin,ymax,N)';
    hx = gx(2)-gx(1);
    hy = gy(2)-gy(1);
    [x,y] = meshgrid(gx,gy);
    bx = MSx(x,y);
    by = MSy(x,y);
    % boundary conditions
    rad = 0.3;
    ia = find((x+1).^2 +y.^2 < rad^2);
    ib = find((x-1).^2 +y.^2 < rad^2);
    t = linspace(0,2*pi,200)';
    cA = [-1 + rad*cos(t),rad*sin(t)];
    cB = [1 + rad*cos(t),rad*sin(t)];
    % Operator L^*: -(div b)mu -bx*Dxmu - by*Dymu + 0.5*epsilon*Delta mu
    % The generator L: bx*Dxmu + by*Dymu + 0.5*epsilon*Delta mu
    N2 = N^2;
    divb = MSdivb(x,y);
    divb = divb(:);
    % set up differentiation matrices for Neumann BCs
    I = speye(N);
    e = ones(N,1);
    E = spdiags([e,-e],[1,-1],N,N);
    E(1,2) = 0;
    E(N,N-1) = 0;
    Dxmatr = 0.5*kron(E,I)/hx;
    Dymatr = 0.5*kron(I,E)/hy;
    E2 = spdiags([e,-2*e,e],-1:1,N,N);
    E2(1,2) = 2; 
    E2(N,N-1) = 2;
    Dxx = kron(E2,I)/hx^2;
    Dyy = kron(I,E2)/hy^2;
    Delta = Dxx+Dyy;
    Bx = spdiags(bx(:),0,N2,N2);
    By = spdiags(by(:),0,N2,N2);
    % find the invariant distribution
%     DivB = spdiags(divb,0,N2,N2);
%     Lstar = -DivB - Bx*Dxmatr - By*Dymatr + 0.5*epsilon*Delta;
    L = Bx*Dxmatr + By*Dymatr + 0.5*epsilon*Delta;
    % use [V,D] = eigs(A,k,'smallestabs')
    [V,eval] = eigs(L',1,'smallestabs');
    fprintf('eval = %d\n',eval);
    summu = sum(V);
    mu = V/(summu*hx*hy);
    mu = reshape(mu,[N,N]);
    if fig == 'y'
        figure;
        hold on
        imagesc(gx,gy,mu);
        colorbar
        set(gca,'Fontsize',fsz);
        plot(cA(:,1),cA(:,2),'color','w','Linewidth',3);
        plot(cB(:,1),cB(:,2),'color','w','Linewidth',3);
        grid
        axis tight
        daspect([1,1,1])
        colorbar
    end
    % the generator for the forward process
    L = Bx*Dxmatr + By*Dymatr + 0.5*epsilon*Delta;
%     % the generator for the time-reversed process
%     aux_x = epsilon*spdiags(Dxmatr*mu(:)./mu(:),0,N2,N2);
%     aux_y = epsilon*spdiags(Dymatr*mu(:)./mu(:),0,N2,N2);
%     ind0 = find(mu < 1e-12);
%     aux_x(ind0) = 0;
%     aux_y(ind0) = 0;
%     LR = -Bx*Dxmatr - By*Dymatr + aux_x*Dxmatr + aux_y*Dymatr + 0.5*epsilon*Delta;

    % compute the forward committor
    q = zeros(N2,1);
    q(ib) = 1;
    Lab = L;
    Lab(ia,:) = 0;
    Lab(ib,:) = 0;
    rhs = -Lab*q;
    rhs(ib) = 1;
    nin = (1:N2)';
    nin([ia;ib]) = [];
    Iab = speye(N2);
    Iab(nin,:) = 0;
    Lab = Lab + Iab;
    Lab(nin,ia) = 0;
    Lab(nin,ib) = 0;

    qplus = reshape(Lab\rhs,[N,N]);
    if fig == 'y'
        figure;
        hold on;
        imagesc(gx,gy,qplus);
        plot(cA(:,1),cA(:,2),'color','w','Linewidth',3);
        plot(cB(:,1),cB(:,2),'color','w','Linewidth',3);
        set(gca,'Fontsize',fsz);
        axis tight
        title('Forward committor','Fontsize',fsz);
        colorbar
        daspect([1,1,1])
    end
    % compute the backward committor
    P = spdiags(mu(:), 0, N2, N2);
    Pinv = spdiags(1./mu(:), 0, N2, N2);
    Lhat = Pinv*L'*P;

    q = zeros(N2,1);
    q(ia) = 1;
    Lab = Lhat;
    Lab(ia,:) = 0;
    Lab(ib,:) = 0;
    rhs = -Lab*q;
    rhs(ia) = 1;
    nin = (1:N2)';
    nin([ia;ib]) = [];
    Iab = speye(N2);
    Iab(nin,:) = 0;
    Lab = Lab + Iab;
    Lab(nin,ia) = 0;
    Lab(nin,ib) = 0;
    qminus = reshape(Lab\rhs,[N,N]);
    if fig == 'y'
        figure;
        hold on
        imagesc(gx,gy,qminus);
        plot(cA(:,1),cA(:,2),'color','w','Linewidth',3);
        plot(cB(:,1),cB(:,2),'color','w','Linewidth',3);
        set(gca,'Fontsize',fsz);
        title('Backward committor','Fontsize',fsz);
        axis tight
        daspect([1,1,1])
        colorbar;
    end
    % compute the reactive current
%     Jx1 = bx.*mu.*qminus-0.5*epsilon*reshape(Dxmatr*(mu(:).*qminus(:)),[N,N]);
%     Jy1 = by.*mu.*qminus-0.5*epsilon*reshape(Dymatr*(mu(:).*qminus(:)),[N,N]);
%     figure;
%     hold on;
%     imagesc(gx,gy,sqrt(Jx1.^2+Jy1.^2));
%     plot(cA(:,1),cA(:,2),'color','w','Linewidth',3);
%     plot(cB(:,1),cB(:,2),'color','w','Linewidth',3);
%     title('Current1','Fontsize',fsz);
%     axis tight
%     colorbar;
    % alternative formula
    Jx = bx.*mu-epsilon*0.5*reshape(Dxmatr*mu(:),[N,N]);
    Jy = by.*mu-epsilon*0.5*reshape(Dymatr*mu(:),[N,N]);
    Jx2 = qplus.*qminus.*Jx +...
        0.5*epsilon*mu.*(qminus.*reshape(Dxmatr*qplus(:),[N,N]) - ...
        qplus.*reshape(Dxmatr*qminus(:),[N,N]));
    Jy2 = qplus.*qminus.*Jy +...
        0.5*epsilon*mu.*(qminus.*reshape(Dymatr*qplus(:),[N,N]) - ...
        qplus.*reshape(Dymatr*qminus(:),[N,N]));
    if fig == 'y'
        figure;
        hold on;
        J = sqrt(Jx2.^2+Jy2.^2);
        Iplot = 4:4:N-3;
        imagesc(gx(Iplot),gy(Iplot),J(Iplot,Iplot));
        plot(cA(:,1),cA(:,2),'color','w','Linewidth',3);
        plot(cB(:,1),cB(:,2),'color','w','Linewidth',3);
        set(gca,'Fontsize',fsz);
        if a == 10
            plot(map(1,:),map(2,:),'Linewidth',3,'color','r');
            plot(map(1,:),-map(2,:),'Linewidth',3,'color','r');
            plot([0,1],[0,0],'Linewidth',3,'color','r');
        else
            if a < 4
                plot([-1,1],[0,0],'Linewidth',3,'color','r');
            end
        end
    %     title('Reactive current','Fontsize',fsz);
        axis tight
        daspect([1,1,1])
        colorbar;
    end
    % transition rate
%    r1_ab = hy*sum(Jx1(:,(N+1)/2));
    r_ab = hy*sum(Jx2(:,(N+1)/2));
    fprintf('r_ab = %d\n',r_ab);
end
%% compute the prefactor
if flag == 1
    N = 2049;
    fname = sprintf("beta%d/solutionAccepted%d.csv",a,N);
    % columns of data:
    % 1 2 3   4      5          6        7       8  9  10 11 12 13   14  15 16    
    % x,y,U,err(U),err(DUx),err(DUy),update_type,Ux,Uy,x1,y1,x2,y2,lambda,a0,a1.
    % (x,y) are the coordinates of the mesh point
    % (x1,y1), (x2,y2) are the parents of (x,y)
    % (xlam,ymal) = (x1,y1) + lambda*((x2,y2) - (x1,y1))
    % a0, a1 = tangent of the angle between the line [(xlam,ylam),(x,y)] and the
    % MAP at the endpoints
    data = readmatrix(fname);
    xmin = -2; xmax = 0;
    ymin = -1; ymax = 1;
    gx = linspace(xmin,xmax,N)';
    gy = linspace(ymin,ymax,N)';
    [x,y] = meshgrid(gx,gy);
    hx = gx(2)-gx(1);
    hy = gy(2)-gy(1);
    h = hx;
    bx = MSx(x,y);
    by = MSy(x,y);
    if abs(h-hy)>1e-15
        fprintf('hx = %d is not equal to hy = %d\n',hx,hy);
        return;
    end
    % read data
     U = reshape(data(1:end-1,3),[N,N]);
%     figure;
%     hold on;
%     imagesc(U);
    
    Ux = reshape(data(1:end-1,7),[N,N]);
    Uy = reshape(data(1:end-1,8),[N,N]);
    px1 = data(1:end-1,10)-1;
    py1 = data(1:end-1,11);
    px2 = data(1:end-1,12)-1;
    py2 = data(1:end-1,13);
    lam = data(1:end-1,14);
    a0 = data(1:end-1,15);
    a1 = data(1:end-1,16);
    % find indices of Accepted points
    ind = find(isfinite(U));
    % equation for the prefactor
    % \nabla * (b + 0.5*\nabla U)C + (b + \nabla U)*\nabla C = 0
    %% compute fac1 := \nabla * (b + 0.5*\nabla U)
    lx = bx + 0.5*Ux;
    ly = by + 0.5*Uy;
    % test rotational component: it must align with level sets of U
    labs = sqrt(lx.^2+ly.^2);
%     figure;
%     hold on
%     contour(gx,gy,U,linspace(0,max(max(U)),30));
%     quiver(x(ind),y(ind),lx(ind)./labs(ind),ly(ind)./labs(ind),'color','k');
%     daspect([1,1,1])
%     axis tight

%    test function Dx
    % bxx = 1-3*x.^2-a*y.^2;
    % bxx1 = Dx(bx,hx);
    % max(max(abs(bxx(:,2:end-1)-bxx1(:,2:end-1))))
    divl = Dx(lx,h) + Dy(ly,h); % divergence of the rotational component
    divlb = divl./sqrt(bx.^2+by.^2);
    if a < 4
        map = [linspace(-1,0,401);zeros(1,401)];
    end
    divlbmap = interp2(x,y,divlb,map(1,:),map(2,:));
    dx = zeros(size(map));
    dx = map - circshift(map,[0,1]); 
    dx(:,1) = zeros(size(map,1),1);
    lxy = cumsum(sqrt(sum(dx.^2,1)));
    l = lxy(end);
    dl = l/(length(x)-1);
    BRint = 2*dl*sum(divlbmap(2:2:end-1)); % midpoint rule
    BRfac = exp(BRint);
    detHeq = 16; % the hessian of the quasipotential at the equilibroum is [4,0;0,4]
    ixstar = N-3;
    if a > 4
        iystar = (N+1)/2+2;
    else
        iystar = (N+1)/2;
    end
    Hxx = 0.5*(Ux(iystar,ixstar+1)-Ux(iystar,ixstar-1))/h
    Hyy = 0.5*(Uy(iystar+1,ixstar)-Uy(iystar-1,ixstar))/h
    Hxy = 0.5*(Uy(iystar,ixstar+1)-Uy(iystar,ixstar-1))/h
    detHstar = abs(Hxx*Hyy-Hxy^2);
    lambdastar = 1; % the positive eval of the Jacobian of b at (0,0)
    ExpTau = (2*pi/lambdastar)*sqrt(detHstar/detHeq)*BRfac*exp(U((N+1)/2,N)/epsilon);
    BRrate = 1/ExpTau;
    fprintf('Bouchet-Reygner rate = %d\n',BRrate);
    pref = (2*pi/lambdastar)*sqrt(detHstar/detHeq)*BRfac;
    Ubar = U((N+1)/2,N);
    r_ab = BRrate;
    fname = sprintf('BRdata_a%d.mat',a);
    save(fname,'pref','Ubar','BRfac');
end
end

%%

function dx = Dx(f,hx)
dx = zeros(size(f));
dx(:,2:end-1) = 0.5*(f(:,3:end)-f(:,1:end-2))/hx;
end

function dy = Dy(f,hy)
dy = zeros(size(f));
dy(2:end-1,:) = 0.5*(f(3:end,:)-f(1:end-2,:))/hy;
end

function plotcubic(p1,p2,zlam,h,zhat,e1hat,e2hat,a0,a1,midpoint,j)
figure(j);
hold on;
grid;
poly = @(t)a0*t-a0*t.^2/h+(a0+a1)*t.^2.*(h-t)/h^2;
t = linspace(0,h,20);
pts = zlam+e1hat*t + e2hat*poly(t);
plot(pts(1,:),pts(2,:));
plot([p1(1),p2(1)],[p1(2),p2(2)]);
plot(midpoint(1),midpoint(2),'.','Markersize',20);
plot(zhat(1),zhat(2),'.','Markersize',20,'color','k');
%daspect([1,1,1]);
end

