function gmam()
makemovie = 0; % if makemovie == 1, make move
% Heymann's and V.E.'s geometric minimum action method
% compites the minimum action path from 
% [-1,0] (asymp. stable equilibrium) to [0,0] (saddle)
% for the Maier & Stein model
% b1 = x-x.^3-10*x.*y.^2;
% b2 = -(1+x.^2).*y;
xi = [-1; 0];
xf = [0; 0];
n = 401;
n1 = n-1;
n2 = n - 2;
%% define the initial path
t = linspace(0,1,n);
x = [xi(1) + (xf(1)-xi(1))*t;sin(pi*t)];
[x,l] = reparametrization(x); % Note: |x'| = l;

%% Set up figure
figure(1);
clf; hold on;
[xx,yy] = meshgrid(-1.5:0.1:1.5,-1:0.1:1);
grid;
% plot vector field
[bb1,bb2] = bfield([xx(:)';yy(:)']);
bb1 = reshape(bb1,size(xx));
bb2 = reshape(bb2,size(xx));
bb = sqrt(bb1.^2 + bb2.^2);
quiver(xx,yy,bb1./(bb+1e-6),bb2./(bb+1e-6));
% plot initial path
mypath = plot(x(1,:),x(2,:),'r','Linewidth',2);
set(gca,'Fontsize',20);
xlabel('x','Fontsize',20);
ylabel('y','fontsize',20);
% setup movie file
if makemovie == 1
    fps = 24; % frames per second
    writerObj = VideoWriter('GMAMdemo.mp4','MPEG-4'); % Name  of movie file
    writerObj.FrameRate = fps; % How many frames per second.
    open(writerObj); 
end
%% Parameters of the GMAM
tau = 1.0e-2; % time step
tol = 1.0e-7; 
%% start
k = 0;
h = (1/n1); % the path is parametrized from 0 to 1
r = tau*n1*n1*ones(n,1);
D2 = spdiags([r -2*r r],-1:1,n2,n2);
nor = Inf;
while nor > tol 
    xold = x;
    dxa = 0.5*(circshift(x,[0,-1])-circshift(x,[0,1]))/h;   % x' along the path
    [b1,b2] = bfield(x);
    lam = sqrt(b1.^2 + b2.^2)/l; % |b|/|x'|
    dlam = 0.5*(circshift(lam,[0,-1])-circshift(lam,[0,1]))/h;
    [b1x,b1y,b2x,b2y] = bgrad(x);
    dbx = (b2x-b1y).*dxa(2,:); % [[(grad b)^T - grad b ]*x'] -- x-component
    dby = (b1y-b2x).*dxa(1,:); % [[(grad b)^T - grad b ]*x'] - y-component
    bx = b1x.*b1 + b2x.*b2; % (grad b)^T b -- x-component
    by = b1y.*b1 + b2y.*b2; % (grad b)^T b -- y-component
    % linearly implicit scheme
    mymatr=(eye(n2) - diag(lam(2:n1).^2)*D2);
    
    rhsx = x(1,2:n1) + tau*(lam(2:n1).*dbx(2:n1) - bx(2:n1) + lam(2:n1).*dlam(2:n1).*dxa(1,2:n1));
    rhsy = x(2,2:n1) + tau*(lam(2:n1).*dby(2:n1) - by(2:n1) + lam(2:n1).*dlam(2:n1).*dxa(2,2:n1));

    rhsx(1) = rhsx(1)+tau*x(1,1)*(n1*lam(2))^2;
    rhsy(1) = rhsy(1)+tau*x(2,1)*(n1*lam(2))^2;
    rhsx(n2) = rhsx(n2)+tau*x(1,end)*(n1*lam(n1))^2;
    rhsy(n2) = rhsy(n2)+tau*x(2,end)*(n1*lam(n1))^2;

    x(1,2:n1) = (mymatr\rhsx')';
    x(2,2:n1) = (mymatr\rhsy')';

    [x,l] = reparametrization(x);
    k = k + 1;
    nor = norm(x - xold)/tau;
    fprintf('iter # %d:  res = %d\n',k,nor);

    figure(1);
    set(mypath,'Xdata',x(1,:),'Ydata',x(2,:));
    drawnow;
    if makemovie == 1
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    end
end
save('MSmap.mat','x');
% [~,isort] = sort(sqrt(sum(x.^2,1)),'ascend');
% [~,Y] = ode45(@maier_stein,[0,1000],x(:,isort(1))');
% plot(Y(:,1),Y(:,2),'b','Linewidth',2);
% [~,Y] = ode45(@maier_stein,[0,1000],x(:,isort(2))');
% plot(Y(:,1),Y(:,2),'b','Linewidth',2);

if makemovie == 1
    close(writerObj); % Saves the movie.
end
end
%%
function [x,l] = reparametrization(x)
% x in an d by n array. Row i of x is the coordinate x_i along the path
% returns a uniformly reparametrized path and its length l
t = linspace(0,1,size(x,2));
dx = zeros(size(x));
dx = x - circshift(x,[0,1]); 
dx(:,1) = zeros(size(x,1),1);
lxy = cumsum(sqrt(sum(dx.^2,1)));
l = lxy(end);
x = interp1(lxy/l,x',t)';
end
%%
function [b1, b2] = bfield(x)
% Input is an d by n array 
b1 = x(1,:) - x(1,:).^3 - 10*x(1,:).*x(2,:).^2;
b2 = -(ones(1,size(x,2)) + x(1,:).^2).*x(2,:);
end
%%
function [b1x, b1y, b2x, b2y] = bgrad(x)
% Input is an d by n array 
b1x = ones(1,size(x,2)) - 3*x(1,:).^2 - 10*x(2,:).^2;
b1y = -20*x(1,:).*x(2,:);
b2x = -2*x(1,:).*x(2,:);
b2y = -(ones(1,size(x,2)) + x(1,:).^2);
end
%%
function dy = maier_stein(~,y)
dy = zeros(2,1);
dy(1) = y(1) - y(1)^3 - 10*y(1)^2*y(2);
dy(2) = -(1 + y(1)^2)*y(2);
end


