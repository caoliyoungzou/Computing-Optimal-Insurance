gam = 0.5; %risk aversion parameter of the insured
lam = 1; % 1/lam is the mean of the loss distribution
f = @(x) lam*exp(-lam*x); % pdf of the loss distribution
%setting up the parameters for the premium principle
alp1 = 0.1;
alp2 = 0.2;
ded1 = 1;
ded2 = 2;


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
integrand2 = @(xx,z) exp(gam*z)*lam*exp(-lam*xx);
integrand3 = @(xx,z) exp(gam*(xx-ded1))*lam*exp(-lam*xx);
integrand4 = @(xx,z) exp(gam*(z+1/gam*log(1+alp1)))*lam*exp(-lam*xx);
integrand5 = @(xx,z) exp(gam*(xx-ded2))*lam*exp(-lam*xx);
integrand6 = @(xx,z) exp(gam*(z+1/gam*log(1+alp1+alp2)))*lam*exp(-lam*xx);

syms z
dfun = @(z) int(integrand1,0,z)+ int(integrand2,xx,z,ded1+z)+int(integrand3,xx,ded1+z,ded1+z+log(1+alp1)/gam)+int(integrand4,xx,ded1+z+log(1+alp1)/gam,ded2+z+log(1+alp1)/gam)+int(integrand5,xx,ded2+z+log(1+alp1)/gam,ded2+z+log(1+alp1+alp2)/gam)+int(integrand6,xx,ded2+z+log(1+alp1+alp2)/gam,100)-exp(gam*z);
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
        %Solution to (2.12)
        elseif x(k) > d(n+1) && x(k) < (ded1+d(n+1))
            I(n+1,k) = x(k)-1/gam*log(M(n));
            sum2 = sum2 + exp(gam*(x(k)-I(n+1,k)))*f(x(k))*dx;
        elseif x(k)>= (ded1+d(n+1)) && x(k)<= (ded1+log(M(n)*(1+alp1))/gam)
            I(n+1,k) = ded1;
            sum2 = sum2 + exp(gam*(x(k)-I(n+1,k)))*f(x(k))*dx;
        elseif x(k)>(ded1+log(M(n)*(1+alp1))/gam) &&  x(k)<=(ded2+log(M(n)*(1+alp1))/gam)
            I(n+1,k) = x(k) -1/gam*log(M(n)*(1+alp1));
            sum2 = sum2 + exp(gam*(x(k)-I(n+1,k)))*f(x(k))*dx;
        elseif x(k)>(ded2+log(M(n)*(1+alp1))/gam) && x(k)<(ded2+log(M(n)*(1+alp1+alp2))/gam)
            I(n+1,k) = ded2;
            sum2 = sum2 + exp(gam*(x(k)-I(n+1,k)))*f(x(k))*dx;
        elseif x(k)>(ded2+log(M(n)*(1+alp1+alp2))/gam)
            I(n+1,k) = x(k) -1/gam*log(M(n)*(1+alp1+alp2));
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

%Producing Figure 3
f3 = figure;
subplot(1,2,1)
nvec = 1:1:n;
plot(nvec,M(1:n),"Linewidth",1.6,'Color','b');
hold on
plot([0,nvec], ones(1,n+1)*M(n+1),"Linewidth",1.6,'Color','k')
ylim([1,1.25])
xlim([0,n])
ax = gca;
ax.FontSize = 14;
xlabel('$n$',"interpreter","latex",'FontSize',20);
ylabel('$M_n$','FontSize',20,"interpreter","latex");
set(gca,'YTick',[1 1.05 1.1 1.15 1.2 1.2288])
subplot(1,2,2)
plot(x,I(n+1,:),"Linewidth",1.6,'Color','b')
ax = gca;
ax.FontSize = 14;
ylim([0,4.5])
xlim([0,5])
xlabel('$x$',"interpreter","latex",'FontSize',20);
ylabel('$\hat{I}_g$',"interpreter","latex",'FontSize',20);
set(gcf,'Position',[300 300 1500 450])