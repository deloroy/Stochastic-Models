// question 1

function [X, T]=simul_clients(mu, lambda, t)
  // simule le comportement de la file jusqu'à l'instant t passé
  // pour mu et lambda les parametres des lois exponentielles des temps de service et d'arrivée :
  // stocke dans T les différents instants de changement d'état (longueur) de la file 
  // stocke dans X les états successifs de la file : X(i) est la longueur de la file pour t=[T(i),T(i+1)[
  X = [0];
  T = [0];
  while T($) < t
    if X($) == 0 then
      T = [T, T($) + grand(1, 1, 'exp', 1/lambda)];
      X = [X, 1];
    else
      T = [T, T($) + grand(1, 1, 'exp', 1/(mu + lambda))];
      X = [X, X($) - 1 + 2 * (rand(1) < (lambda / (mu + lambda)))];
    end
  end
endfunction

//Simulations du comportement de la chaine pour différentes densités de trafic
t = 100;

lambda = 1; // paramètre pour la loi du temps d'arrivée
mu = 1; // paramètre pour la loi du temps de service
[X, T]=simul_clients(mu, lambda, t);
clf();plot2d2(T, X);xtitle('Une trajectoire de la longueur de file. Rho=1.');halt();

lambda = 1;
mu = 2;
[X, T]=simul_clients(mu, lambda, t);
clf();plot2d2(T, X);xtitle('Une trajectoire de la longueur de file. Rho=1/2.');halt();

lambda = 2;
mu = 1;
[X, T]=simul_clients(mu, lambda, t);
clf();plot2d2(T, X);xtitle('Une trajectoire de la longueur de file. Rho=2.');halt();

// ---------------------------------------------------------------------------
// question 2

function [E_er, V_er,T] = moments_ergodique(mu, lambda, t)
    // simule le comportement de la file jusqu'à l'instant t passé
    // et retourne les espérances et les variances de file calculées le long de la trajectoire jusqu'à chaque changement d'état inférieur à t
    // en les vecteurs E_er et V_er
    // E_er(i) est l'estimation de l'espérance calculée jusqu'au temps t=T(i+1); 
  [X, T] = simul_clients(mu, lambda, t);
  T_diff = T(2:$) - T(1:$-1);
  E_er=[0];
  V_er=[0];
  for t=1:length(T_diff)
    E_er = [E_er,E_er($)+X(t)*T_diff(t)];
    V_er = [V_er,sum((X(1:t) - E_er($)/T(t+1)) ** 2 .* T_diff(1:t))]; 
  end
  E_er(2:$) = E_er(2:$)./T(2:$);
  V_er(2:$) = V_er(2:$)./T(2:$);
endfunction

lambda = 1;
mu = 2;
rho = lambda / mu;
t = 1000;

// valeurs théoriques de l'espérance et de la variance de la longueur de file en régime stationnaire
E_th = rho / (1-rho);
V_th = rho / (1-rho) ** 2;

// valeurs obtenues en faisant l'intégrale de ces quantités sur une trajectoire
[E_er, V_er,T] = moments_ergodique(mu, lambda, t);

// vérification du théorème ergodique en régime stationnaire
E_er($) - E_th
V_er($) - V_th

// affichage des courbes des espérances et variances le long de la trajectoire en fonction du nombre d'itérations
clf();plot2d2(T(2:$), E_er(2:$));xtitle('Espérance calculée le long de trajectoire arrêté au temps t. Rho=1/2.');halt();
clf();plot2d2(T(2:$), V_er(2:$));xtitle('Variance calculée le long de trajectoire arrêté au temps t. Rho=1/2.');halt();

// courbe des ecarts finaux E_er($) - E_th et V_er($) - V_th en fonction de rho
vect_rho = 0.01:0.01:0.99;
ecarts_E = []; //ecart final E_er($) - E_th fonction de rho
ecarts_V = []; //ecart final V_er($) - V_th fonction de rho
for k=1:length(vect_rho)
    rho=vect_rho(k);
    E_th = rho / (1-rho);
    V_th = rho / (1-rho) ** 2;
    [E_er, V_er,T] = moments_ergodique(1, rho, t);
    ecarts_E = [ecarts_E, abs(E_er($)-E_th)];
    ecarts_V = [ecarts_V, abs(V_er($)-V_th)];
end
clf();plot2d2(vect_rho, ecarts_E);xtitle('Erreur espérance théorique - estimée ergodique au bout de 1000 itérations, en fonction de la densité de trafic rho');halt();
clf();plot2d2(vect_rho, ecarts_V);xtitle('Erreur variance théorique - estimée ergodique au bout de 1000 itérations, en fonction de la densité de trafic rho');halt();
    
    

// ---------------------------------------------------------------------------
// question 3

function [dist_er] = distribution_ergodique(mu, lambda, t)
    // simule le comportement de la file jusqu'à l'instant t passé
    // et retourne les distributions de la longueur de file calculées le long de la trajectoire arrêtée à différents instants jusqu'à t
    // dist_er(k,i) est l'estimation de probabilité que la longueur de la file soit égale à i-1, calculée le long de la trajectoire arrêtée jusqu'à t/k
  [X, T] = simul_clients(mu, lambda, t);
  T_diff = T(2:$) - T(1:$-1);
  X_max = max(X);
  dist_er = zeros(10, X_max+1); //10 estimations de la distribution au fil de la trajectoire jusqu'à t
  dt = floor(length(T)/10.); 
  for k=1:10
    for i = 1:X_max+1
       dist_er(k, i) = sum(bool2s(X(1:k*dt-1) == i-1) .* T_diff(1:k*dt-1)) / T(k*dt);
    end
  end
endfunction

lambda = 1;
mu = 2;
rho = lambda / mu;

// distribution de la longueur de file obtenue en faisant l'intégrale des probabilités sur une trajectoire
dist_er = distribution_ergodique(mu, lambda, t);

// distribution théorique de la longueur de file en régime stationnaire
dist_th = (rho**(0:(size(dist_er,2) - 1)))*(1-rho);

// affichage des deux distributions (vérification du théorème ergodique)
dist = zeros(size(dist_er,1)+1,size(dist_er,2));
dist(1:size(dist_er,1),:) = dist_er;
dist($,:) = dist_th;

clf();
xtitle('Distributions de la longueur de file calculées le long de trajectoire arrêtée à 10 instants espacés jusqu à T_final, puis distribution théorique (dernières barres). Rho=1/2.');
bar(0:(size(dist,2) - 1), dist');
halt();

// ---------------------------------------------------------------------------
// question subsidiaire

function [X, T]=simul_clients_multiguich(mu, lambda, t,K)
  // pour K guichets
  // simule le comportement de la file jusqu'à l'instant t passé
  // pour mu et lambda les parametres des lois exponentielles des temps de service et d'arrivée :
  // stocke dans T les différents instants de changement d'état (longueur) de la file 
  // stocke dans X les états successifs de la file : X(i) est la longueur de la file pour t=[T(i),T(i+1)[
  X = [0];
  T = [0];
  while T($) < t
    if X($) == 0 then
      T = [T, T($) + grand(1, 1, 'exp', 1/lambda)];
      X = [X, 1];
    else
      T = [T, T($) + grand(1, 1, 'exp', 1/(min(X($),K)*mu + lambda))];
      X = [X, X($) - 1 + 2 * (rand(1) < (lambda / (min(X($),K)*mu + lambda)))];
    end
  end
endfunction

//Simulations du comportement de la chaine pour différentes densités de trafic
t = 100;
K = 10;

lambda = 5; // paramètre pour la loi du temps d'arrivée
mu = 1; // paramètre pour la loi du temps de service
[X, T]=simul_clients_multiguich(mu, lambda, t, K);
clf();plot2d2(T, X);xtitle('Une trajectoire de la longueur de file pour 10 guichets. Rho=1/2.');halt();

lambda = 20;
mu = 1;
[X, T]=simul_clients_multiguich(mu, lambda, t, K);
clf();plot2d2(T, X);xtitle('Une trajectoire de la longueur de file pour 10 guichets. Rho=2.');halt();

lambda = 10;
mu = 1;
[X, T]=simul_clients_multiguich(mu, lambda, t, K);
clf();plot2d2(T, X);xtitle('Une trajectoire de la longueur de file pour 10 guichets. Rho=1.');halt();
