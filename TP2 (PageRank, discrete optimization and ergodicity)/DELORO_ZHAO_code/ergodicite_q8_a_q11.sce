// --------------------------------------------------------------------------
// Question 8
// Test de "vérification" du théorème ergodique : 
// L'espérance du cout calculée avec le PageRank comme la limite du cout moyen
// sur une trajectoire de page en page infiniment longue

// D�finition de la fonction �crite en question 5 :
// Calul de pi par convergence de la suite p_{k+1}=P'*p_{k}
function [pi]=pi_iterative(n, P)  
  p=ones(n,1);  
  while %t  
    pn = P' * p; 
    if norm(pn-p,%inf) < 10*%eps then break;end
    p = pn;
  end  
  pi= pn' / sum(pn); 
endfunction  

// Fonction coût 
function y=r(x)  
  y=x^2  
endfunction  
 
n=4;  
P=rand(n,n)  
pr=sum(P,'c');  
P = P ./ (pr*ones(1,n));  


function cerg=ergodique_markov_T(T,P)
  // Calcule le coût moyen d'une page le long d'une trajectoire à T pas
  // Choix aléatoire de la page de départ : loi uniforme sur les pages
  idx_i = int(grand(1, 1, "unf", 1, 5));
  cerg = r(idx_i);
  for t = 1:T-1
    // Choix aléatoire de la page suivante
    // en suivant la loi de probabilité donnée par P(idx_i, :)
    proba_j = cumsum(P(idx_i, :));
    idx_j = rand(1, 1, 'uniform');
    [m, idx_i] = max(bool2s(idx_j < proba_j));
    cerg = cerg + r(idx_i); 
  end
  cerg = cerg / T;
endfunction  
 
function [cerg,pi]=ergodique_markov(P)
  // Calcule l'espérance du cout d'une page avec le PageRank 
  // (limite du théorème ergodique)
  pi = pi_iterative();
  cerg = 0;
  for i = 1:n
    cerg = cerg + pi(1, i) * r(i)
  end
endfunction  
 
// Test (vérification du théorème ergodique)
T=100000; CT=ergodique_markov_T(T,P);  
[c,pi]=ergodique_markov(P);  
c-CT

// --------------------------------------------------------------------------
// Question 9
// Vérification que l'espérance du cout calculée avec le PageRank s'obtient
// en résolvant un certain système linéaire

// Le noyau de P-I est engendré par ones(n,1)  
[x0,K]=linsolve(P- eye(n,n),zeros(n,1));

// le projecteur spectral sur Espace propre associé a 1  
Pr = ones(n,1)*pi; // [pi;pi;pi;....]  
A = P-eye(n,n);    // A -Id  
S = Pr - inv(Pr-A) // Pr-A est inversible  
// vérifier que S*Pr et Pr*S sont nuls  
clean(S*Pr)  
clean(Pr*S)  
// A*w  + R - c= 0  
// A*c         = 0  
R=r((1:n)');  
// vérifions que w=-S*R et c=Pr*R sont solution du systeme linaire  
w= -S*R;  
c= Pr*R;  
A*w + R -c  
A*c  
// Noter que w n'est pas unique, on peut rajouter à w les elts du noyau de A  
// Montrons inversement que c doit être egal à Pr*R  
// Pr*A est nul  
Pr*A  
// on doit donc avoir  
// Pr*R - Pr*c = 0 et A*c =0  
// en sommant  
// Pr*R = (Pr-A)*c  
// c = (Pr-A)^-1 *Pr*R  
// c = (Pr-S)*Pr*R = Pr*Pr*R -S*Pr*R = Pr*R  
// car Pr est un projecteur Pr^2 = Pr et S*Pr = 0  
clean(Pr^2-Pr)  
clean(S*Pr)  
// conclusion c doit valoir Pr*R  
// on le vérifie avec linsolve  
 
[x0,K]=linsolve([A,-eye(n,n);zeros(n,n),A],[R;zeros(n,1)]);  
// on vérifie bien que e = Pr*R

// --------------------------------------------------------------------------
// Cas particulier : P a la forme de la matrice de Google

P1=rand(n,n)  
pr=sum(P1,'c');  
P1 = P1 ./ (pr*ones(1,n));  
 
z=grand(1,n,'unf',0,1);  
z=z/sum(z);  
 
alpha = 0.8;  
 
P = alpha*P1 + (1-alpha)*ones(n,1)*z;  
 
// Les couts Rm(i,j)  
Rm = grand(n,n,'unf',0,1);

// --------------------------------------------------------------------------
// Question 10
// Calcul de w par résolution d'un système linéaire

// On le vérifie umeriquement  
// trouver la solution de  
// w = alpha*P1*w + sum(P.*Rm,'c')  

w = linsolve(alpha * P1 - eye(n, n), sum(P.*Rm, 'c'));
 
// calcul de c  
c = (1-alpha)*z*w  
 
// (w,c) solution du pb ergodique ?  
w + c - (P*w + sum(P.*Rm, 'c'))  
 
// --------------------------------------------------------------------------
// Question 11
// Calcul de w par convergence de suite

function w=iterative_c(tol)
  b = sum(P.*Rm, 'c');
  w = ones(n, 1);  
  while %t  
    wn = alpha * P1 * w + b;
    if norm(wn-w,%inf) < 10*%eps then break;end
    w = wn;
  end  
  w = wn
endfunction  
 
w=iterative_c(10*%eps);  
// calcul de c  
c = (1-alpha)*z*w  
 
// (w,c) solution du pb ergodique ?  
w + c - (P*w + sum(P.*Rm,'c'))

