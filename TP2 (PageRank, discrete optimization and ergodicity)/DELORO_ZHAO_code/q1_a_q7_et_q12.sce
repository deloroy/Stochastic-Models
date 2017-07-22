n=10; // Nombre de pages  
 
function show_adj(Adj,diameters)  
  [lhs,rhs]=argn(0);  
  if rhs < 2 then diameters = 30*ones(1,n);end  
  graph = mat_2_graph(sparse(Adj),1,'node-node');  
  graph('node_x')=300*cos(2*%pi*(1:n)/(n+1));  
  graph('node_y')=300*sin(2*%pi*(1:n)/(n+1));  
  graph('node_name')=string([1:n]);  
  graph('node_diam')= diameters;  
  // graph('node_color')= 1:n;  
  // show_graph(graph);  
  rep=[1 1 1 1 2 2 2 2 2 2 2 2 2];  
  plot_graph(graph,rep);  
endfunction  
 
Adj=grand(n,n,'bin',1,0.2);show_adj(Adj);halt();

// --------------------------------------------------------------------------
// Question 2
// Construction de la matrice de transition P  
// associÃ©e Ã  une matrice d'adjacence.  
// Pss: transition d'origine,  
// P: matrice de google  
// z: vecteur de teleportation  
// d: vecteur vaut 1 si le degrÃ© vaut zero et 0 sinon  

function [P,Pss,Pprim,d,z,alpha]=google(Adj)

  z = rand(1, n);
  z = z / sum(z);

  alpha = 0.85
  
  P = zeros(n, n);
  Pss = zeros(n, n);
  Pprim = zeros(n, n); 
  for i = 1:n
    if sum(Adj(i, :)) ~= 0
      Pss(i, :) = Adj(i, :) / sum(Adj(i, :));
      Pprim(i, :) = Pss(i, :);
    else
      Pprim(i, :) = z;
    end
  end

  P = alpha * Pprim + (1-alpha) * ones(n, 1) * z;
  d = bool2s(sum(Adj, 2) == 0);

endfunction  
 
[P,Pss,Pprim,d,z,alpha]=google(Adj); 
 
// verification que P est stochastique   
sum(P,'c')

// --------------------------------------------------------------------------
// Question 3 
// Calcul de P'*x en conservant le caractÃ¨re creux des opÃ©rations
x=rand(n,1);  
y1=P'*x;
y2=alpha*Pss'*x + z' * ones(1, n) * x + alpha * z' * (d' - 1) * x;
y1-y2 // VÃ©rification du calcul 1.0D-15  

// --------------------------------------------------------------------------
// Question 4
// Calul de pi (Page Rank des pages du graphe) par Ã©tude du spectre de P'
[U, lambda] = spec(P');  
// U matrice des vecteurs propres, lambda matrice diagonale (termes : valeurs propres)
idx = 0;
// Recherche de la valeur propre 1 dans lambda
for i=1:n do
  value_diff = abs(lambda(i, i) - 1);
  if value_diff < 1D-5 then
    idx = i;
    break
  end
end
// Vecteur propre associÃ© Ã  la valeur propre 1
pi = real(U(:, idx)'); 
pi = pi / sum(pi);

disp("The value of pi calculated with spectrum of P");
disp(pi);
xbasc();
// Dessin du graphe en choisissant le diamÃ¨tre des cercles des noeuds avec pi
show_adj(Adj,int(300*pi)); 
halt();

// --------------------------------------------------------------------------
// Question 5
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
 
pi=pi_iterative(n, P);  
clean(pi*P - pi); // VÃ©rification que pi est mesure invariante associÃ©e Ã  la chaine
disp("The value of pi calculated by convergence of sequence");
disp(pi);

// --------------------------------------------------------------------------
// Question 6
// Calul de pi par convergence de la suite p_{k+1}=P'*p_{k} 
// En utilisant des matrices sparses
function [pi]=pi_iterative_sparse(n, Pss, d, z, alpha)  
  p=ones(n,1);  
  while %t  
    pn = alpha*Pss'*p + z' * ones(1, n) * p + alpha * z' * (d' - 1) * p; 
    if norm(pn-p,%inf) < 10*%eps then break;end
    p = pn;
  end   
  pi= pn' / sum(pn); 
endfunction  
 
pi=pi_iterative_sparse(n, Pss, d, z, alpha);  
clean(pi*P - pi) // VÃ©rification
disp("The value of pi calculated by sparse iteration");
disp(pi);

// --------------------------------------------------------------------------
// Question 7
// RÃ©solution du problÃ¨me d'optimisation du PageRank par optimisation discrÃ¨te

function [P,Pss,Pprim,d]=google2(Adj,z, alpha)
    // Dans l'objectif de pouvoir comparer la matrice d'adjacence optimisée à la matrice initiale :
    // redéfinition d'une fonction "google", qui crée les mêmes matrices P,Pss,Pprim et d
    // A partir du vecteur de téléportation z et du alpha qui sont cette fois passés en entrée
    // (autrement dit la matrice de transition P est connue avant "google2", pas d'aléa dans cette fonction)
  P = zeros(n, n);
  Pss = zeros(n, n);
  for i = 1:n
    if sum(Adj(i, :)) ~= 0
      Pss(i, :) = Adj(i, :) / sum(Adj(i, :));
      Pprim(i, :) = Pss(i, :);
    else
      Pprim(i, :) = z;
    end
  end

  P = alpha * Pprim + (1-alpha) * ones(n, 1) * z;
  d = bool2s(sum(Adj, 2) == 0);
endfunction

m = n/2;
p = 2;

function [Adj_optim, cost]=pagerank(Adj,z,alpha, m, p, n)
    //Retourne la matrice d'adjacence optimisant le PageRank des m-premieres pages
    //Ainsi que la valeur du PageRank des m-premières pages en la variable cost
  Adj_optim = Adj;
  [P,Pss,Pprim,d]=google2(Adj_optim, z, alpha);
  pi=pi_iterative_sparse(n, Pss, d, z, alpha);
  cost = sum(pi(1:m));
  
  // Nombre total de combinaisons Ã  tester :
  // PrÃ©sence ? d'un lien de chacune des p premiÃ¨res pages Ã  chacune des pages m+1...n
  test_num = 2 ** (p * (n-m)); 
  
  for i=1:test_num do 
    idx_bin = dec2bin(i-1); // Une combinaison possible
    
    // Construction de la matrice d'adjacence y correspondant
    Adj_0 = Adj;
    Adj_0(1:p, m+1:$) = zeros(p, n-m);
    for j=1:length(idx_bin)
      idx_x = 1 + floor(j / (n-m+1)); //de quel page part-on ?
      idx_j = m + modulo(j, n-m+1); //vers quel page va-t'on ?
      if part(idx_bin, j) == '0' then
        Adj_0(idx_x, idx_j) = 0;
      else
        Adj_0(idx_x, idx_j) = 1;
      end
    end
    
    // Calcul du PageRank
    [P,Pss,Pprim,d]=google2(Adj_0, z, alpha);
    pi=pi_iterative_sparse(n, Pss, d, z, alpha); // PageRank (pi)
    cost_0 = sum(pi(1:m)); // PageRank des m premiÃ¨res pages

    // Obtient-on un meilleur PageRank des m premiÃ¨res pages ?
    if cost_0 > cost then
        Adj_optim = Adj_0;
        cost = cost_0;
    end
  end
endfunction

[Adj_optim, cost]=pagerank(Adj, z,alpha,m, p, n);
[P_optim,Pss_optim,Pprim_optim,d_optim]=google2(Adj_optim, z, alpha);
pi_optim=pi_iterative_sparse(n, Pss_optim, d_optim, z, alpha);

// Affichage du graphe "optimisé"
show_adj(Adj_optim,int(300*pi_optim)); 

// --------------------------------------------------------------------------
// Questions 8 à 11 : voir l'autre fichier de code

// --------------------------------------------------------------------------
// Question 12
// Maximisation du PageRank via ergodicitÃ©

Rm1 = ones(m, n);
Rm0 = zeros(n-m, n);
Rm = [Rm1; Rm0];

function [Adj_optim, cost]=pagerank_ergodic(Adj, z, alpha, m, p, n)
  w = ones(n, 1);
  w_before = zeros(n,1);
  Adj_optim = Adj;
  [P_optim,Pss_optim,Pprim_optim,d_optim]=google2(Adj_optim, z, alpha);
   
  while norm(w - w_before, %inf) > %eps
      // sauvegarde du précédent w
      w_before=w;
      
      //1. Choix du meilleur contrôle Adj_optim
      
      cost = alpha * Pprim_optim * w + sum(P_optim.*Rm,2);  
      
      for k=1:p
        // Choix du contrôle sur les liens allant du noeud k vers les noeuds m+1,...,n
        // On veut maximiser la k-ème composante de cost
        cost_k = cost(k,1);
        
        // Nombre total de combinaisons Ã  tester :
        // PrÃ©sence ? d'un lien de la i-ème page Ã  chacune des pages m+1...n
        test_num = 2 ** (n-m); 
      
        for i=1:test_num do 
          idx_bin = dec2bin(i-1); // Une combinaison possible
          
          // Construction de la matrice d'adjacence y correspondant
          Adj_0 = Adj;
          Adj_0(k, m+1:$) = zeros(1, n-m);
          for j=1:length(idx_bin) //vers quel page va-t'on ?
            if part(idx_bin, j) == '0' then
              Adj_0(k, m+j) = 0;
            else
              Adj_0(k, m+j) = 1;
            end
          end
        
          // Calcul de la k-ième composante du cout à optimiser
          [P_0,Pss_0,Pprim_0,d_0]=google2(Adj_0, z, alpha);
          cost_0 = alpha * Pprim_0 * w + sum(P_0.*Rm,2);
          cost_0_k = cost_0(k,1);

          // Obtient-on un meilleur PageRank des m premiÃ¨res pages ?
          if cost_0_k > cost_k then
              Adj_optim(k,:) = Adj_0(k,:);
              cost_k=cost_0_k;
          end
        end
     
      //disp(Adj_optim);
   
      end 
    
      //2. Calcul du cout ergodique associé à ce meilleur controle via w$
       
      [P_optim,Pss_optim,Pprim_optim,d_optim]=google2(Adj_optim, z, alpha);
      [w] = linsolve([alpha * Pprim_optim - eye(n,n)] , [sum(P_optim.*Rm,2)]);
       
   end
    
   cost = (1-alpha)*z*w;
   
endfunction

[Adj_optim_2, cost_2]=pagerank_ergodic(Adj, z,alpha,m, p, n);
[P_optim_2,Pss_optim_2,Pprim_optim_2,d_optim_2]=google2(Adj_optim_2, z, alpha);
pi_optim_2=pi_iterative_sparse(n, Pss_optim_2, d_optim_2, z, alpha);

// Affichage du graphe "optimisé"
show_adj(Adj_optim_2,int(300*pi_optim_2)); 

