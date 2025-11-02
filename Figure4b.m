tic
%Parameters
a=75; 
b=10^5;
r=0.023;
m=28740;
mu1=21.05; 
mu2=0.347;
mu3=0.2;
mu4=0.72;
p1=2.3;
p2=0.0125;
p3=(3.7)*10^(-6);
p4=1.44*10^(-5);
d0=1.032*10^4;
gamma=9.12;
beta= 5.2;
eta=7.2*10^(-2);
k=10^9;

for j=1:1:100
Initial=[ 0 1*j*10^6 d0/mu2 (gamma*d0/mu2)/(mu3+(p4*eta*d0)/(mu2*mu4)) eta*d0/(mu2*mu4)]; 
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);


tspan1 = [0,0+2/24];
for i=1:6
m=28740;
[t,y] = ode23s(@(t,y) odefcn(t,y,a,b,r,m,mu1,mu2,mu3,mu4,p1,p2,p3,p4,d0,beta,gamma,eta,k), tspan1, y0, opt);
semilogy(t,y(:,2),'k')
s=semilogy(t,y(:,2),'k');
set(s,'linewidth',2.5);
xlabel('Time (Days)');
ylabel('# Cells');
hold on

y0 = y(end,:);
tspan1 = [7*(i-1)+2/24,7*i];
m=0;

[t,y] = ode23s(@(t,y) odefcn(t,y,a,b,r,m,mu1,mu2,mu3,mu4,p1,p2,p3,p4,d0,beta,gamma,eta,k), tspan1, y0, opt);
semilogy(t,y(:,2),'k')
s=semilogy(t,y(:,2),'k');
set(s,'linewidth',2.5);


hold on
y0 = y(end,:);
tspan1 = [7*i,7*i+2/24];
end
 
hold on 
y0 = y(end,:);

for i=7:75
    
    if mod(i+1,4)==0 && i<32 && i>10
       m=28740;
[t,y] = ode23s(@(t,y) odefcn(t,y,a,b,r,m,mu1,mu2,mu3,mu4,p1,p2,p3,p4,d0,beta,gamma,eta,k), tspan1, y0, opt);
semilogy(t,y(:,2),'k')
s=semilogy(t,y(:,2),'k');
set(s,'linewidth',2.5);
hold on


y0 = y(end,:);
tspan1 = [7*(i-1)+2/24,7*i];
m=0;


[t,y] = ode23s(@(t,y) odefcn(t,y,a,b,r,m,mu1,mu2,mu3,mu4,p1,p2,p3,p4,d0,beta,gamma,eta,k), tspan1, y0, opt);
semilogy(t,y(:,2),'k')
s=semilogy(t,y(:,2),'k');
set(s,'linewidth',2.5);


hold on
y0 = y(end,:);
tspan1 = [7*i,7*i+2/24];
    else
m=0;
[t,y] = ode23s(@(t,y) odefcn(t,y,a,b,r,m,mu1,mu2,mu3,mu4,p1,p2,p3,p4,d0,beta,gamma,eta,k), tspan1, y0, opt);
semilogy(t,y(:,2),'k')
s=semilogy(t,y(:,2),'k');
set(s,'linewidth',2.5);
hold on

y0 = y(end,:);
tspan1 = [7*(i-1)+2/24,7*i];
m=0;

[t,y] = ode23s(@(t,y) odefcn(t,y,a,b,r,m,mu1,mu2,mu3,mu4,p1,p2,p3,p4,d0,beta,gamma,eta,k), tspan1, y0, opt);
semilogy(t,y(:,2),'k')
s=semilogy(t,y(:,2),'k');
set(s,'linewidth',2.5);


hold on
y0 = y(end,:);
tspan1 = [7*i,7*i+2/24];


    end
end
end

hold on
%yline(10^9, 'color', [0.5 0.5 0.5],'LineWidth', 3, 'linestyle', '--' ,'Label', 'Carrying capacity of tumor');
yline(10^9, 'color', [0.8500, 0.3250, 0.0980]		, 'LineWidth', 5, 'linestyle', '--', ...
    'Label', 'Carrying capacity of tumor', ...
    'LabelHorizontalAlignment', 'right');

legend('Tumor cells (T)')
fontsize(20,"points")
fontweight='bold';
xlabel('Days','FontWeight', 'bold')
ylabel('# Cells','FontWeight', 'bold')
ax = gca;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
yticks([ 1e-1 1e0 1e2 1e4 1e6 1e8 1e10 ])               % specify positions
yticklabels({'10^{-1}','10^0','10^2','10^4','10^6','10^8','10^{10}'})  % specify labels
%ax.YAxis.Exponent = 2;
xlim([0 500])
ylim([0.1 1*10^10])
xticks(0:50:500)
set(gca,'YminorTick','off')
xlabel('Time (Days)','FontWeight', 'bold')
ylabel('# Cells','FontWeight', 'bold')
hold off
%-------------------------------------------
toc




function dydt =odefcn(t,y,a,b,r,m,mu1,mu2,mu3,mu4,p1,p2,p3,p4,d0,beta,gamma,eta,k)
dydt = zeros(5,1);
M=y(1);
T=y(2);
D=y(3);
E=y(4);
R=y(5);

dydt = [ -mu1*M+m;
r*(1-M/(M+a))*T*(1-T/k)-p1*T*M/(M+a)-p3*E*T*exp((-R*T/k)*(1-M/(M+a))); 
d0-mu2*D+beta*(1-D/b)*T*p1*M/(M+a); 
gamma*D-mu3*E-p4*E*R;
eta*D*(1-M/(M+a))-mu4*R-p2*R*M/(M+a)];
end



