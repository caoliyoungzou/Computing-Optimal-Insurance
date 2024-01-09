lam = 1; %1/lam is the mean of the loss distribution
f = @(x) lam*exp(-lam*x); %pdf of the loss distribution
gam = 2; %risk aversion parameter of the insured
alp = 0.5; %premium loading


N = 100;
K = 50000;
x = linspace(0,15,K+1);
dx = x(2)-x(1);
d = zeros(1,N+1); %vector to store d_n
M = zeros(1,N+1); %vector to store M_n
I = zeros(N,K+1); %vector to store I_n

%Calculating the theoretic value of M
syms xx
integrand1 = @(xx) exp(gam*xx)*lam*exp(-lam*xx);
Ifun = @(xx,z)(1/gam*lambertw(gam/2/alp*exp(gam*(xx-z+1/2/alp)))-1/2/alp);
integrand2 = @(xx,z) exp(gam*(xx-Ifun(xx,z)))*lam*exp(-lam*xx);
syms z
dfun = @(z) int(integrand1,0,z)+ int(integrand2,xx,z,100)-exp(gam*z);
dnum = eval(vpasolve(dfun(z) == 0,z,1));
M_theo = exp(gam*dnum);

%Applying Algorithm 2.1
%Step 1. Initialize M_0
M(1) = 1; 
for n = 1:N
    %Step 2. Compute d_n
    d(n+1) = 1/gam*log(M(n));
    sum1 = 0; %First term in eqn (2.15)
    sum2 = 0; %Second term in eqn (2.15)
    %Step 3. Compute I_n(x)
    for k = 1:(K+1)
        if x(k) <= d(n+1) 
            I(n+1,k) = 0;
            sum1 = sum1 + exp(gam*(x(k)))*f(x(k))*dx;
        else 
            I(n+1,k) = 1/gam*lambertw(gam/2/alp*exp(gam*(x(k)-d(n+1)+1/2/alp)))-1/2/alp; %Solution to eqn (2.12)
            sum2 = sum2 + exp(gam*(x(k)-I(n+1,k)))*f(x(k))*dx;
        end
    end
    %Step 4. Compute M_n
    M(n+1) = sum1 + sum2;
    %Step 5. Convergence of M
    if abs(M(n+1)-M(n)) < 10^-6
        break
    end
end

%Producing Figure 2
f2 = figure;
subplot(1,2,1)
nvec = 1:1:n;
plot(nvec,M(1:n),"Linewidth",1.6,'Color','b');
hold on
plot([0,nvec], ones(1,n+1)*M(n+1),"Linewidth",1.6,'Color','k')
ylim([1,5.5])
xlim([0,n])
ax = gca;
ax.FontSize = 14;
xlabel('$n$',"interpreter","latex",'FontSize',20);
ylabel('$M_n$','FontSize',20,"interpreter","latex");
set(gca,'YTick',[1 1.5 2 2.5 3 3.5 4 4.5 5 5.4214])
subplot(1,2,2)
plot(x(:),I(n+1,:),"Linewidth",1.6,'Color','b');
ax = gca;
ax.FontSize = 14;
xlabel('$x$',"interpreter","latex",'FontSize',20);
ylabel('$\hat{I}_g$',"interpreter","latex",'FontSize',20);
ylim([0,3.5])
xlim([0,5])
set(gcf,'Position',[300 300 1500 450])


