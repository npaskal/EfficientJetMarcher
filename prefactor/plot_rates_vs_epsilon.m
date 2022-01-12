function plot_rates_vs_epsilon()
fsz = 24; % fontsize
ep = (0.01:0.01:0.1)'; % values of epsilon: dX = b(X)dt +\sqrt{\epsilon}dW
Neps = length(ep);
% Set the parameter a for the Maier-Stein SDE
a = 3;
if a == 3
    J = 2:Neps;
else
    J = 1:Neps;
end
% Load or generate TPT rates
fname = sprintf('TPTrate_a%d.mat',a);
if ~isfile(fname)
    rate_ab = zeros(Neps,1);
    flag = 0;
    for j = 1 : Neps 
        rate_ab(j) = transition_rate(flag,ep(j),a);
    end
    save(fname,'rate_ab','ep');
else
    TPTdata = load(fname);
    rate_ab = TPTdata.rate_ab;
end
% Load data file with Bouchet-Reygner prefactor and quasipotential barrier
fname = sprintf('BRdata_a%d.mat',a);
if ~isfile(fname)
    r = transition_rate(1,ep(end),a);
end
BRdata = load(fname);
pref = BRdata.pref;
Ubar = BRdata.Ubar;
BRtau = pref*exp(Ubar./ep);
BRfac = BRdata.BRfac;
fprintf('BRfac = %d\n',BRfac);
% Plot expected exit time vs 1/epsilon
figure;
hold on;
plot(1./ep(J),1./rate_ab(J),'Linewidth',1,'Marker','s','Markersize',6,'Displayname','TPT');
plot(1./ep(J),BRtau(J),'Linewidth',1,'Marker','<','Markersize',6,'Displayname','Bouchet-Reyner');
set(gca,'Fontsize',fsz,'YScale','log');
xlabel('1/\epsilon','Fontsize',fsz);
ylabel('E[\tau_{AB}]','Fontsize',20);
legend
% Plot prefactors vs epsilon
figure;
hold on;
plot(ep(J),exp(-Ubar./ep(J))./rate_ab(J),'Linewidth',1,'Marker','s','Markersize',6,'Displayname','TPT');
plot(ep(J),pref*ones(size(ep(J))),'Linewidth',1,'Marker','<','Markersize',6,'Displayname','Bouchet-Reyner');
set(gca,'Fontsize',fsz,'YScale','log');
xlabel('\epsilon','Fontsize',fsz);
ylabel('E[\tau_{AB}]exp(-U(0,0)/\epsilon)','Fontsize',20);
save('rate_data_a10.mat','rate_ab','BRtau','ep');
legend
end

