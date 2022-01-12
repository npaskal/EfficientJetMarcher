function MSprefactor()
% if flag == 0, compute the prefactor and find pdf = Cexp(-U/epsilon)
% if flag == 1, compute pdf by solving forward Kolmogorov equation
close all
a = 3; % the parameter a in the Maier-Stein SDE
epsilon = 0.1; % the value of epsilon in dX = b(X)dt + \sqrt{\epsilon}dW
%%
fsz = 24;
cfig = 5;
cvals = 0.2:0.4:10;
%% Maier-Stein vector field
MSx = @(x,y)x-x.^3-a*x.*y.^2;
MSy = @(x,y)-(1+x.^2).*y;
Pot = @(x,y)x.^4/4-x.^2/2+x.^2.*y.^2/2+y.^2/2;
MSdivb = @(x,y)-4*x.^2-a*y.^2;
for flag = 0:1
%% compute the prefactor
if flag == 0
   fname = sprintf('WKB_eps%.1f_a%d.mat',epsilon,a);
   if isfile(fname)
        % save('WKB_eps0.1_a3.mat','gx','gy','WKBmu','C','U','epsilon','N')
        WKBdata = load(fname);
        N = WKBdata.N;
        gx = WKBdata.gx;
        gy = WKBdata.gy;
        [x,y] = meshgrid(gx,gy);
        hx = gx(2)-gx(1);
        hy = gy(2)-gy(1);
        h = hx;
        WKBmu = WKBdata.WKBmu;
        C = WKBdata.C;
        U = WKBdata.U;
    else
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
        figure;
        hold on
        contour(gx,gy,U,linspace(0,max(max(U)),30),'Linewidth',1);
        set(gca,'Fontsize',fsz);
        grid
        axis tight
        xlabel('x','Fontsize',fsz);
        ylabel('y','Fontsize',fsz);
        axis([-1.5,0,-0.75,0.75]);
        daspect([1,1,1])

    %    test function Dx
        % bxx = 1-3*x.^2-a*y.^2;
        % bxx1 = Dx(bx,hx);
        % max(max(abs(bxx(:,2:end-1)-bxx1(:,2:end-1))))
        divl = Dx(lx,h) + Dy(ly,h); % divergence of the rotational component
        % central index
        icenter = (N+1)/2;
        C = zeros(N);
        k_init = round(N/32);
        init = icenter-k_init:icenter+k_init;
        ind_init = find(U<0.001);
        C(ind_init) = 1;
        [Usort,update_order] = sort(U(ind),'ascend');
        mark = zeros(N);
        Naccepted = length(ind);
        mark(ind_init) = 1; % points with computed C

        for j = 1 : Naccepted
            if mod(j,1000) == 0
                fprintf('j = %d\n',j);
            end
            index = ind(update_order(j)); % index of current point
            if mark(index) == 0 % if index is not computed yet
                ip1 = getindex(px1(index),py1(index),h,N);
                ip2 = getindex(px2(index),py2(index),h,N);
                if mark(ip1) == 1 && mark(ip2) == 1               
                    [C,mark] = compute_prefactor(C,mark,data,x,y,divl,MSx,MSy,bx,by,index,ip1,ip2,px1,py1,px2,py2,lam,a0,a1);
                else 
                    fprintf('U(index) = %d, U(ip1) = %d, U(ip2) = %d\n',U(index),U(ip1),U(ip2));
                    if U(ip1) > U(ip2)
                        index = ip1;
                        ip1 = getindex(px1(index),py1(index),h,N);
                        ip2 = getindex(px2(index),py2(index),h,N);
                        [C,mark] = compute_prefactor(C,mark,data,x,y,divl,MSx,MSy,bx,by,index,ip1,ip2,px1,py1,px2,py2,lam,a0,a1);
                    else
                        index = ip2;
                        ip1 = getindex(px1(index),py1(index),h,N);
                        ip2 = getindex(px2(index),py2(index),h,N);
                        [C,mark] = compute_prefactor(C,mark,data,x,y,divl,MSx,MSy,bx,by,index,ip1,ip2,px1,py1,px2,py2,lam,a0,a1);
                    end 
                    index = ind(update_order(j));
                    ip1 = getindex(px1(index),py1(index),h,N);
                    ip2 = getindex(px2(index),py2(index),h,N);
                    [C,mark] = compute_prefactor(C,mark,data,x,y,divl,MSx,MSy,bx,by,index,ip1,ip2,px1,py1,px2,py2,lam,a0,a1);                
                end
            end
        end
        WKBmu = zeros(N);
        WKBmu(ind) = C(ind).*exp(-U(ind)/epsilon);
    end 
    ifinite = find(isfinite(WKBmu));
    Z = hx*hy*sum(WKBmu(ifinite));
    WKBmu = WKBmu/Z;
%     Z = hx*hy*sum(sum(WKBmu));
%     WKBmu = WKBmu/Z;
%     WKBmusum = sum(WKBmu(:))
%     WKBmu = WKBmu/(WKBmusum*h^2);
%     sum(WKBmu(:))*h^2
    kplot = round((N-1)/256);
    Iplot = kplot:kplot:N-kplot+1;
    figure(1);
    hold on
%    imagesc(gx(Iplot),gy(Iplot),WKBmu(Iplot,Iplot));
    contour(x,y,WKBmu,cvals,'Linewidth',1);
    colorbar
    set(gca,'Fontsize',fsz);
    grid
    axis tight 
    xlabel('x','Fontsize',fsz);
    ylabel('y','Fontsize',fsz);
    axis([-1.5,0,-0.75,0.75]);
    daspect([1,1,1])
    
    figure(2);
    hold on
    imagesc(gx(Iplot),gy(Iplot),C(Iplot,Iplot));
    colorbar
    set(gca,'Fontsize',fsz);
    plotMAPs(a);
    grid
    axis tight  
    xlabel('x','Fontsize',fsz);
    ylabel('y','Fontsize',fsz);
    axis([-1.5,0,-0.75,0.75]);
    daspect([1,1,1])
    if a == 10
        caxis([0,2.6]);
    end
    fname = sprintf('WKB_eps%.1f_a%d.mat',epsilon,a);
    save(fname,'gx','gy','WKBmu','C','U','epsilon','N')
end
%% compute invariant pdf by solving (L^*)mu = 0, hx*hy*sum(mu) = 1
if flag == 1
    N = 1025;
    xmin = -2; xmax = 0;
    ymin = -1; ymax = 1;
    gx = linspace(xmin,xmax,N)';
    gy = linspace(ymin,ymax,N)';
    hx = gx(2)-gx(1);
    hy = gy(2)-gy(1);
    [x,y] = meshgrid(gx,gy);
    bx = MSx(x,y);
    by = MSy(x,y);
    % Operator L^*: -(div b)mu -bx*Dxmu - by*Dymu + 0.5*epsilon*Delta mu
    % The generator L: bx*Dxmu + by*Dymu + 0.5*epsilon*Delta mu
    N2 = N*N;
    bx = bx(:);
    by = by(:);
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
    Bx = spdiags(bx,0,N2,N2);
    By = spdiags(by,0,N2,N2);
    L = Bx*Dxmatr + By*Dymatr + 0.5*epsilon*Delta;
    % use [V,D] = eigs(A,k,'smallestabs')
    [V,eval] = eigs(L',1,'smallestabs');
    fprintf('eval = %d\n',eval);
    summu = sum(V);
    mu = V/(summu*hx*hy);
    mu1 = reshape(mu,[N,N]);
    figure(1);
    Iplot = 4:4:N-3;
    hold on
    % imagesc(gx(Iplot),gy(Iplot),mu1(Iplot,Iplot));
    contour(gx,gy,mu1,cvals,'--','Linewidth',1);
    colorbar
    set(gca,'Fontsize',fsz);
    grid
    axis tight
    xlabel('x','Fontsize',fsz);
    ylabel('y','Fontsize',fsz);
    axis([-1.4,-0.3,-0.5,0.5]);
    daspect([1,1,1])

    %contour(gx,gy,mu1,linspace(0,max(max(mu)),30));
end
end
end


%%
function ind = getindex(x,y,h,N)
jx = round((x+2)/h);
jy = round((y+1)/h);
ind = jx*N + jy + 1;
%fprintf('x = %d, y = %d, jx = %d, jy = %d, ind = %d\n',x,y,jx,jy,ind);
end

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

%%
function [C,mark] = compute_prefactor(C,mark,data,x,y,divl,MSx,MSy,bx,by,index,ip1,ip2,px1,py1,px2,py2,lam,a0,a1)
   % proceed to calculation of the prefactor
    % find xlambda
    Clam = C(ip1) + lam(index)*(C(ip2)-C(ip1));
    % C(x) = C(xlam)*exp(-\int_0^hlam div(l)*sqrt(1+tan^2)/||b||dh
    % approximate the integral using Simpson's rule
    % restore cubic curve approximating MAP segment
    p1 = [px1(index);py1(index)];
    p2 = [px2(index);py2(index)];
    zlam = p1+lam(index)*(p2-p1);
    zhat = [data(index,1)-1;data(index,2)];
    hlam = norm(zhat-zlam); 
    e1hat = zhat-zlam;
    e1hat = e1hat/norm(e1hat);
    e2hat = [-e1hat(2);e1hat(1)];
    poly = @(t)a0(index)*t-a0(index)*t.^2/hlam+(a0(index)+a1(index))*t.^2.*(hlam-t)/hlam^2;
%                midpoint = 0.5*(zhat+zlam)+e2hat*0.125*hlam*(a0(index)-a1(index));
    midpoint = zlam+e1hat*0.5*hlam + e2hat*poly(0.5*hlam);
    divlS = interp2(x,y,divl,[zlam(1),midpoint(1)],[zlam(2),midpoint(2)]);
    blam = norm(MSx(zlam(1),zlam(2)),MSy(zlam(1),zlam(2)));
    bmid = norm(MSx(midpoint(1),midpoint(2)),MSy(midpoint(1),midpoint(2)));   
%     % Simpson's rule
%     Integral = (hlam/6)*(divlS(1)*sqrt(1+a0(index)^2)/blam +...
%         divlS(2)*sqrt(1+0.0625*(a0(index)+a1(index))^2)/bmid +...
%         divl(index)*sqrt(1+a1(index)^2)/norm([bx(index);by(index)]));
%     % midpoint rule
%     Integral = hlam*divlS(2)*sqrt(1+0.0625*(a0(index)+a1(index))^2)/bmid;
%    C(index) = Clam*exp(-Integral); 
    % backward Euler
    C(index) = Clam/(1+hlam*divl(index)/norm([bx(index);by(index)]));
    mark(index) = 1;
end

%%
function plotMAPs(a)
MSx = @(x,y)x-x.^3-a*x.*y.^2;
MSy = @(x,y)-(1+x.^2).*y;
MSJ11 = @(x,y)1-3*x.^2-a*y.^2;
MSJ12 = @(x,y)-2*a*x*y;
MSJ21 = @(x,y)-2*x.*y;
MSJ22 = @(x,y)-1-x.^2;
eq = [-1;0];
% the Jacobian for system linearized at (-a,-a+a^3/3)
J = [-2, 0;0,-2];
A = [J,eye(2);zeros(2),-J']; % The right-hand side of the linearized Hamiltonian system
% find unstable manifold of the equilibrium for the Hamiltonian system
[V,E] = eig(A);
evals = diag(E);
revals = real(evals);
[rsort,jsort] = sort(revals,'descend');
V = V(:,jsort);
if abs(imag(evals(1))) > 1e-15
    MX = [real(V(1:2,1)),imag(V(1:2,1))];
    MP = [real(V(3:4,1)),imag(V(3:4,1))];
else
    MX = V(1:2,1:2);
    MP = V(3:4,1:2);
end
M = (MX'\MP')';

% M = MP*MX^{-1}
% initial conditions for MAPs
rad = 1e-9;
ni = 50;
t = linspace(0,2*pi,ni+1);
t(end) = [];
% the RHS for the Hamiltonian system
fun = @(t,y)[y(3)+MSx(y(1),y(2));y(4)+MSy(y(1),y(2));...
    -MSJ11(y(1),y(2))*y(3)-MSJ21(y(1),y(2))*y(4);...
    -MSJ12(y(1),y(2))*y(3)-MSJ22(y(1),y(2))*y(4)];
options = odeset("AbsTol",1e-12,"RelTol",1e-12,"Events",@events1);
for j = 1 : ni
    z = rad*[cos(t(j));sin(t(j))];
    p = M*z; % initial conditions for momenta on the unstable manifold for the Hamiltonian system 
    x = [z+eq;p];
    [~,Y] = ode45(fun,[0,1000],x,options);
    plot(Y(:,1),Y(:,2),'color','w');
end
end
function [value,isterminal,direction] = events1(t,y)
value(1) = abs(y(1) + 1) - 1;
value(2) = abs(y(2)) - 1;
isterminal(1) = 1;
isterminal(2) = 1;
direction(1:2) = 0;
end
