%% Rosenbrock Function [Punto II]
f0 =@(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
gradf0 =@(x) [400*x(1)^3-400*x(1)*x(2)+2*x(1)-2;200*(x(2)-x(1)^2)];
Hessf0 =@(x) [1200*x(1)^2-400*x(2)+2,-400*x(1);-400*x(1),200];
x0 = [1.2;1.2];
x1 = [-1.2;1];
opt = ones(2);
rho = 0.5;
c = 10.^(-4);
kmax = 100000;
tolgrad = 10.^(-12);
backtrackmax = 50;
tol_pcg = 10.^(-6);
maxit_pcg = 50000;

% tic ;
% [xsol,fsol,gradfsol,k,rate,normgrad]=mod_newton_backtrack(x0,f0,gradf0,Hessf0,kmax,tolgrad,c,rho,backtrackmax,opt);
% xsol
% k
% gradfsol
% time = toc ;
% time 
% 
% tic ;
% [xsol1,fsol1,gradfsol1,k1,rate1,normgrad1]=tru_newton_backtrack_cg(x0,f0,gradf0,Hessf0,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% xsol1
% k1
% gradfsol1
% time1 = toc ;
% time1
% 
% tic ;
% [xsol2,fsol2,gradfsol2,k2,rate2,normgrad2]=tru_newton_backtrack_pcg(x0,f0,gradf0,Hessf0,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% xsol2
% k2
% gradfsol2
% time2 = toc ;
% time2


%% Chained Rosenbrock Function [Punto III - Primo Problema]
f_ch = @(x) ch_ros(x);
Gradf_ch = @(x) grad_ch_ros(x);
Hessf_ch = @(x) hess_ch_ros(x);
n = 1000 ;
x0_ch = ones(n, 1); % Initialize the vector with all elements as 1
x0_ch(1:2:end) = -1.2; % Set odd indices to -1.2
x1_ch=rand(n,1)+0.5; %Random tra 0.5 e 1.5
rho = 0.5;
c = 10.^(-4);
kmax = 10000;
tolgrad = 10.^(-6);
backtrackmax = 50;
opt = ones(n,1) ;
tol_pcg = 10.^(-6);
maxit_pcg = 50000;

% tic ;
% [xsol_ch,fsol_ch,gradfsol_ch,k_ch,rate_ch,normgrad_ch]=mod_newton_backtrack(x0_ch,f_ch,Gradf_ch,Hessf_ch,kmax,tolgrad,c,rho,backtrackmax,opt);
% xsol_ch
% k_ch
% gradfsol_ch
% time_ch = toc ;
% time_ch 
% 
% tic ;
% [xsol_ch1,fsol_ch1,gradfsol_ch1,k_ch1,rate_ch1,normgrad_ch1]=tru_newton_backtrack_cg(x0_ch,f_ch,Gradf_ch,Hessf_ch,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% xsol_ch1
% k_ch1
% gradfsol_ch1
% time_ch1 = toc ;
% time_ch1
% 
% tic ;
% [xsol_ch2,fsol_ch2,gradfsol_ch2,k_ch2,rate_ch2,normgrad_ch2]=tru_newton_backtrack_pcg(x0_ch,f_ch,Gradf_ch,Hessf_ch,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% xsol_ch2
% k_ch2
% gradfsol_ch2
% time_ch2 = toc ;
% time_ch2


%% Problema 76 [Punto III - Secondo Problema]
f76 = @(x) f_76(x);
Gradf76 = @(x) grad_76(x);
Hessf76 = @(x) hess_76(x);
rho = 0.5;
c = 10.^(-4);
n = 1000;
kmax = 10000;
tolgrad = 10.^(-4);
backtrackmax = 50;
x0_76 = 2*ones(n,1);
rng(123)
x1_76= 4*rand(n,1)-2; %Random tra -0.5 e 0.5
x2_76= rand(n, 1) + 3
tol_pcg = 10.^(-6);
maxit_pcg = 5000;
opt=zeros(n,1);

% tic ;
% [xsol_76,fsol_76,gradfsol_76,k_76,rate_76,normgrad_76]=mod_newton_backtrack(x1_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,opt);
% xsol_76
% k_76
% gradfsol_76
% time_76 = toc ;
% time_76 

% [xsol_76,fsol_76,gradfsol_76,k_76,rate_76,normgrad_76]=mod_newton_backtrack(x0_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,opt);
% normgrad1=gradfsol_76
% rate0=rate_76
% 
% [xsol_76,fsol_76,gradfsol_76,k_76,rate_76,normgrad_76]=mod_newton_backtrack(x1_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,opt);
% normgrad2=gradfsol_76
% rate1=rate_76
% 
% [xsol_76,fsol_76,gradfsol_76,k_76,rate_76,normgrad_76]=mod_newton_backtrack(x2_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,opt);
% normgrad3=gradfsol_76
% rate2=rate_76
% 
% plot(rate0)
% hold on
% plot(rate1)
% plot(rate2)
% xlabel('Iterations')
% ylabel('Rate of convergence')
% legend('x0','x1','x2')
% hold off



tic ;
[xsol_761,fsol_761,gradfsol_761,k_761,rate_761,normgrad_761]=tru_newton_backtrack_cg(x2_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
xsol_761
k_761
gradfsol_761
time_761 = toc ;
time_761

% [xsol_761,fsol_761,gradfsol_761,k_761,rate_761,normgrad_761]=tru_newton_backtrack_cg(x0_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad1=normgrad_761
% rate0=rate_761
% 
% [xsol_761,fsol_761,gradfsol_761,k_761,rate_761,normgrad_761]=tru_newton_backtrack_cg(x1_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad2=normgrad_761
% rate1=rate_761
% 
% [xsol_761,fsol_761,gradfsol_761,k_761,rate_761,normgrad_761]=tru_newton_backtrack_cg(x2_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad3=normgrad_761
% rate2=rate_761
% 
% plot(rate0)
% hold on
% plot(rate1)
% plot(rate2)
% xlabel('Iterations')
% ylabel('Rate of convergence')
% legend('x0','x1','x2')
% hold off
% 
% tic ;
% [xsol_762,fsol_762,gradfsol_762,k_762,rate_762,normgrad_762]=tru_newton_backtrack_pcg(x1_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% k_762
% mean(xsol_762)
% gradfsol_762
% time_762 = toc ;
% time_762

% [xsol_762,fsol_762,gradfsol_762,k_762,rate_762,normgrad_762]=tru_newton_backtrack_pcg(x0_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad1=normgrad_762
% rate0=rate_762
% 
% [xsol_762,fsol_762,gradfsol_762,k_762,rate_762,normgrad_762]=tru_newton_backtrack_pcg(x1_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad2=normgrad_762
% rate1=rate_762
% [xsol_762,fsol_762,gradfsol_762,k_762,rate_762,normgrad_762]=tru_newton_backtrack_pcg(x2_76,f76,Gradf76,Hessf76,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad3=normgrad_762
% rate2=rate_762
% plot(rate0)
% hold on
% plot(rate1)
% plot(rate2)
% xlabel('Iterations')
% ylabel('Rate of convergence')
% legend('x0','x1','x2')
% hold off

%% Problema 82 [Punto III - Terzo Problema]
f82 = @(x) f_82(x);
Gradf82 = @(x) grad_82(x);
Hessf82 = @(x) hess_82(x);
n = 1000;
x0_82 = 1/2*ones(n, 1); % Initialize the vector with all elements as 1/2
rng(123)
x1_82=4*rand(n,1)-2;
x2_82=rand(n,1)+1;
x3_82=rand(n,1)+10;
x4_82=rand(n,1)+20;%Random tra -1 e 1
rho = 0.5;
c = 10.^(-4);
kmax = 10000;
tolgrad = 10.^(-10);
backtrackmax = 50;
opt = zeros(n,1) ;
tol_pcg = 10.^(-6);
maxit_pcg = 5000;

tic ;
[xsol_82,fsol_82,gradfsol_82,k_82,rate_82,normgrad_82]=mod_newton_backtrack(x0_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,opt);
k_82
sol=mean(xsol_82)
gradfsol_82
time_82 = toc;
time_82 

% [xsol_82,fsol_82,gradfsol_82,k_82,rate_82,normgrad_82]=mod_newton_backtrack(x0_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,opt);
% normgrad1=gradfsol_82
% rate0=rate_82
% 
% [xsol_82,fsol_82,gradfsol_82,k_82,rate_82,normgrad_82]=mod_newton_backtrack(x1_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,opt);
% normgrad2=gradfsol_82
% rate1=rate_82
% 
% [xsol_82,fsol_82,gradfsol_82,k_82,rate_82,normgrad_82]=mod_newton_backtrack(x2_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,opt);
% normgrad3=gradfsol_82
% rate2=rate_82

% plot(normgrad1)
% hold on
% plot(normgrad2)
% plot(normgrad3)
% xlabel('Iterations')
% ylabel('Norm of the gradient')
% legend('x0','x1','x2')
% hold off


% plot(rate0)
% hold on
% plot(rate1)
% plot(rate2)
% xlabel('Iterations')
% ylabel('Rate of convergence')
% legend('x0','x1','x2')
% hold off

% tic ;
% [xsol_821,fsol_821,gradfsol_821,k_821,rate_821,normgrad_821]=tru_newton_backtrack_cg(x1_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% k_821
% mean(xsol_821)
% normgrad1 = gradfsol_821
% time_821 = toc ;
% time_821
% [xsol_821,fsol_821,gradfsol_821,k_821,rate_821,normgrad_821]=tru_newton_backtrack_cg(x0_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad1=normgrad_821
% rate0=rate_821
% 
% [xsol_821,fsol_821,gradfsol_821,k_821,rate_821,normgrad_821]=tru_newton_backtrack_cg(x1_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad2=normgrad_821
% rate1=rate_821
% 
% [xsol_821,fsol_821,gradfsol_821,k_821,rate_821,normgrad_821]=tru_newton_backtrack_cg(x3_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad3=normgrad_821
% rate2=rate_821
% 
% plot(normgrad1)
% hold on
% plot(normgrad2)
% plot(normgrad3)
% xlabel('Iterations')
% ylabel('Norm of the gradient')
% legend('x0','x1','x3')
% hold off


% plot(rate0)
% hold on
% plot(rate1)
% plot(rate2)
% xlabel('Iterations')
% ylabel('Rate of convergence')
% legend('x0','x1','x3')
% hold off

% tic ;
% [xsol_822,fsol_822,gradfsol_822,k_822,rate_822,normgrad_822]=tru_newton_backtrack_pcg(x1_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% k_822
% mean(xsol_822)
% gradfsol_822
% time_822 = toc ;
% time_822

% [xsol_822,fsol_822,gradfsol_822,k_822,rate_822,normgrad_822]=tru_newton_backtrack_pcg(x0_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad1=normgrad_822
% rate0=rate_822
% 
% [xsol_822,fsol_822,gradfsol_822,k_822,rate_822,normgrad_822]=tru_newton_backtrack_pcg(x1_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad2=normgrad_822
% rate1=rate_822
% 
% [xsol_822,fsol_822,gradfsol_822,k_822,rate_822,normgrad_822]=tru_newton_backtrack_pcg(x4_82,f82,Gradf82,Hessf82,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxit_pcg,opt);
% normgrad3=normgrad_822
% rate2=rate_822
% 
% plot(normgrad1)
% hold on
% plot(normgrad2)
% plot(normgrad3)
% xlabel('Iterations')
% ylabel('Norm of the gradient')
% legend('x0','x1','x4')
% % hold off
% 
% plot(rate0)
% hold on
% plot(rate1)
% plot(rate2)
% xlabel('Iterations')
% ylabel('Rate of convergence')
% legend('x0','x1','x4')
% hold off

%% Modified Newton [Punto I - d]

function [xsol,fsol,gradfsol,k,rate,normgrad]=mod_newton_backtrack...
    (x0,f,gradf,Hessf,kmax,tolgrad,c,rho,backtrackmax,opt)
xsol = x0; %Inizializzazione
k=1;
error_k_0 = norm(x0-opt,2); %Norma 2 errore
normgrad(1)=(norm(gradf(xsol)));
gradfsol = zeros(20,1)
rate=zeros(20,1)
while k<kmax && norm(gradf(xsol)) > tolgrad %Condizione di stop
    backtrack=0;
    alpha = 1; %Passo
    flag = 1;
    Hessfk = Hessf(xsol); %Hessiana soluzione
    beta = norm(Hessfk,"fro"); %Norma di frobenius
    if min(diag(Hessfk)) > 0 %Calcolo tau
        tau = 0;
    else
        tau = beta/2;
    end
    while flag ~= 0 %Finche Hess+E non è positiva definita
        E = tau*eye(length(x0)); %Calcolo E tau*I
        [~,flag] = chol(Hessfk+E); %Flag=0 se riesco a decomporla
        tau = max(2*tau,beta/2);
    end
    B = Hessfk + E; %Calcolo la B
    L=chol(B);
    p=L\(L'\-gradf(xsol));
    while backtrack < backtrackmax && (f(xsol+alpha*p)>f(xsol)+...
            c*alpha*gradf(xsol)'*p)  %Armijo-Goldstein condition
        alpha = rho*alpha;
        backtrack = backtrack+1;
    end
    error_k_1 = norm(xsol-opt,2); %Norma errore prima di aggiornare xsol
    xsol = xsol + p*alpha; %Calcolo xsol
    fsol=f(xsol); %Calcolo la f
    gradfsol(k)=norm(gradf(xsol)); %Calcolo il gradiente
    k=k+1;
    error_k = norm(xsol-opt,2);  %Calcolo il nuovo errore
    if k>=2
        %Calcolo il rate di convergenza
        rate(k-1)=log(error_k/error_k_1)/log(error_k_1/error_k_0); 
        error_k_0 = error_k_1;
    end
    %Salvo i valore della norma del gradiente
    normgrad(k)=(norm(gradf(xsol))); 
end
end


%% Truncated Newton Senza Precondizionamento [Punto I - d]

function [xsol,fsol,gradfsol,k,rate,normgrad]=tru_newton_backtrack_cg(x0,f,gradf,Hessf,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxiterations_pcg,opt)
xsol = x0; %Inizializzazione
k=1;
error_k_0 = norm(x0-opt,2); %Norma 2 errore
normgrad= zeros(50,1)
rate=zeros(50,1)
normgrad(1)=(norm(gradf(xsol)));
while k<kmax && norm(gradf(xsol)) > tolgrad %Condizione di stop
    alpha = 1; %Passo
    backtrack=0;
    Hessfk = Hessf(xsol); %Hessiana soluzione
    Gradfk = gradf(xsol); %gradiente soluzione
    nu=min(0.5,norm(Gradfk)^1/2);
    z=zeros(length(x0),1);
    r=Gradfk;
    d=-r;
    w=0;
    while w < maxiterations_pcg && norm(r)/norm(Gradfk) > tol_pcg
        if d'*Hessfk*d<=0 %Controllo che Hess sia definita positiva o meno
            if w==0 %Se non è definita positiva al primo passo tengo -grad come direzione ed esco
                p=-Gradfk;
                break
            else
                p=z;
                break
            end
        end
        v=Hessfk*d;
        teta=r'*r/(d'*v); %Continuo il calcolo del coniugate gradient
        z=z+teta*d;
        r1=r+teta*v;
        if norm(r1)<norm(Gradfk)*nu %Controllo se è rispettata la norma
            p=z;
            break
        end
        beta=r1'*r1/(r'*r);
        d=-r1+beta*d;
        r=r1;
        w=w+1;
    end
    while backtrack < backtrackmax && (f(xsol+alpha*p)>f(xsol)+c*alpha*Gradfk'*p)  %Armijo-Goldstein condition
        alpha = rho*alpha;
        backtrack = backtrack+1;
    end
    error_k_1 = norm(xsol-opt,2); %Norma errore prima di aggiornare xsol
    xsol = xsol + p*alpha; %Calcolo xsol
    fsol=f(xsol); %Calcolo la f
    gradfsol=norm(gradf(xsol)); %Calcolo il gradiente
    k=k+1;
    error_k = norm(xsol-opt,2);  %Calcolo il nuovo errore
    if k>=2
        rate(k-1)=log(error_k/error_k_1)/log(error_k_1/error_k_0); %Calcolo il rate di convergenza
        error_k_0 = error_k_1;
    end
    normgrad(k)=(norm(gradf(xsol)));
end
end


%% Truncated Newton Con Precondizionamento [Punto I - d]

function [xsol,fsol,gradfsol,k,rate,normgrad]=tru_newton_backtrack_pcg(x0,f,gradf,Hessf,kmax,tolgrad,c,rho,backtrackmax,tol_pcg,maxiterations_pcg,opt)
xsol = x0; %Inizializzazione
k=1;
error_k_0 = norm(x0-opt,2); %Norma 2 errore
normgrad=zeros(150,1)
rate=zeros(150,1)
normgrad(1)=(norm(gradf(xsol)));
while k<kmax && norm(gradf(xsol)) > tolgrad %Condizione di stop
    alpha = 1; %Passo
    backtrack=0;
    Hessfk = Hessf(xsol); %Hessiana soluzione
    Gradfk = gradf(xsol); %gradiente soluzione
    nu=min(0.1,norm(Gradfk)^1/2);
    z=zeros(length(x0),1);
    [L,U]=ilu(sparse(Hessfk));
    r=Gradfk;
    y=U\(L\r);
    d=-y;
    w=0;
    while w < maxiterations_pcg && norm(r)/norm(Gradfk) > tol_pcg
        if d'*Hessfk*d<=0 %Controllo che Hess sia definita positiva o meno
            if w==0 %Se non è definita positiva al primo passo tengo -grad come direzione ed esco
                p=-Gradfk;
                break
            else
                p=z;
                break
            end
        end
        v=Hessfk*d;
        teta=r'*y/(d'*v); %Continuo il calcolo del coniugate gradient
        z=z+teta*d;
        r1=r+teta*v;
        if norm(r1)<norm(Gradfk)*nu %Controllo se è rispettata la norma
            p=z;
            break
        end
        y1=U\(L\r1);
        beta=r1'*y1/(r'*y);
        d=-y1+beta*d;
        r=r1;
        y=y1;
        w=w+1;
        p = z;
    end
    while backtrack < backtrackmax && (f(xsol+alpha*p)>f(xsol)+c*alpha*Gradfk'*p)  %Armijo-Goldstein condition
        alpha = rho*alpha;
        backtrack = backtrack+1;
    end
    error_k_1 = norm(xsol-opt,2); %Norma errore prima di aggiornare xsol
    xsol = xsol + p*alpha; %Calcolo xsol
    fsol=f(xsol); %Calcolo la f
    gradfsol=norm(gradf(xsol)); %Calcolo il gradiente
    k=k+1;
    error_k = norm(xsol-opt,2);  %Calcolo il nuovo errore
    if k>=2
        rate(k-1)=log(error_k/error_k_1)/log(error_k_1/error_k_0); %Calcolo il rate di convergenza
        error_k_0 = error_k_1;
    end
    normgrad(k)=(norm(gradf(xsol))); %Salvo i valore della norma del gradiente
end
end


%% Definizione Chained Rosenbrock Function
function y = ch_ros(x)
    n = length(x); % Dimensione Punti
    y=0;
    for i = 2:n
        y= y + 100*(x(i-1)^2-x(i))^2 + (x(i-1)-1)^2;
    end
end

function grad = grad_ch_ros(x)
    n = length(x);
    grad = zeros(n,1);
    grad(1) = 400*x(1)^3-400*x(1)*x(2)+2*x(1)-2;
    for i = 2:n-1
        grad(i) = -200*x(i-1)^2+200*x(i)+400*x(i)^3-400*x(i)*x(i+1)+2*x(i)-2;
    end
    grad(n) = 200*x(n) - 200*x(n-1)^2;
end

function hess = hess_ch_ros(x)
    n = length(x);
    hess = zeros(n,n);
    hess(1,1) = 1200*x(1)-400*x(2)+2;
    hess(1,2) = -400*x(1);
    hess(n,n) = 200;
    hess(n,n-1) = hess(n-1,n);
    for i=2:n-1
        hess(i,i) = 200+1200*x(i)^2-400*x(i+1)+2;
        hess(i,i+1) = -400*x(i);
        hess(i,i-1) = -400*x(i-1);
    end
end


%% Definizione Funzione Problema 76
function y = f_76(x)
    n = length(x); % Dimensione Punti
    y=1/2*(x(1)-x(2)^2/10)^2;
    for i=2:n-1
        y= y+1/2*(x(i)-x(i+1)^2/10)^2;
    end
    y=y+1/2*(x(n)-x(1)^2/10)^2;
end

function grad = grad_76(x)
    n = length(x);
    grad = zeros(n,1);
    grad(1) = 1/50*x(1)^3-x(2)^2/10-x(1)*x(n)/5+x(1);
    for i = 2:n-1
        grad(i) = 1/50*x(i)^3-x(i+1)^2/10-x(i)*x(i-1)/5+x(i);
    end
    grad(n) = 1/50*x(n)^3-x(1)^2/10-x(n)*x(n-1)/5+x(n);
end

function hess = hess_76(x)
    n = length(x);
    hess = zeros(n,n);
    hess(1,1) = 1+3/50*x(1)^2-x(n)/5;
    hess(1,2) = -x(2)*50;
    hess(1,n) = -x(1)/5;
    hess(n,1) = -x(1)/5;
    hess(n,n) = 3/50*x(n)^2-x(n-1)/5+1;
    hess(n,n-1) = -x(n)/5;
    for i=2:n-1
        hess(i,i) = 3/50*x(i)^2 - x(i-1)/5+1 ;
        hess(i,i+1) = -x(i+1)/5;
        hess(i,i-1) = -x(i)/5;
    end
end


%% Definizione Problema 82
function y = f_82(x)
    n = length(x); % Dimensione Punti
    y=1/2*x(1)^2;
    for i = 2:n
        y=y+1/2*(x(i)+cos(x(i-1))-1)^2;
    end
end

function grad = grad_82(x)
    n = length(x);
    grad = zeros(n,1);
    grad(1)=x(1)+sin(x(1))-sin(x(1))*x(2)-sin(x(1))*cos(x(1));
    for i = 2:n-1
        grad(i)=-1+x(i)+cos(x(i-1))+sin(x(i))-sin(x(i))*x(i+1)-sin(x(i))*cos(x(i));
    end
    grad(n)=-1+x(n)+cos(x(n-1));
end

function hess = hess_82(x)
    n = length(x);
    hess=zeros(n,n);
    hess(1,1)=1+cos(x(1))-cos(x(1))*x(2)-cos(x(1))^2+sin(x(1))^2;
    hess(1,2)=-sin(x(1));
    hess(n,n)=1;
    hess(n,n-1)=-sin(x(n-1));
    for i=2:n-1
        hess(i,i) = 1+cos(x(i))-cos(x(i))*x(i+1)-cos(x(i))^2+sin(x(i))^2;
        hess(i,i+1) = -sin(x(i));
        hess(i,i-1) = -sin(x(i-1));
    end
end