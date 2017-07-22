// Notations
// A : matrice de transition de S (taille 2x2)
// B : proba conditionnelle de Y sachant S (taille 2x4)
// filt : proba conditionnelle de S sachant Y (taille 4x2)
// base : nombre entre 1 et 4 indexant la base (ACGT)

// Lit le fichier ADN (NE PAS MODIFIER)
function [adn]=lecture_adn(toto)
  A = mgetl(toto);
  A=[A;'tccggtgatc cgacaggtta cg'];
  B = ascii(A);
  B(find(B==32))=[]; 
  // replace ACGT by 1234
  code= "acgt" ;  
  acode = ascii(code);
  val = [1,2,3,4]; 
  for i=1:length(code);
  B(find(B==acode(i))) = val(i);
  end;
  adn=B;
endfunction;

// Programme principal (NE PAS MODIFIER)
// Retourne la probabilité d'être en chaque valeur possible de l'état adjoint
// connaissant la séquence d'ADN, ainsi que
// l'évolution des matrices A et B au cours des itérations
function [region, p0i, Aevol,Bevol]=result(adn,nbre_etat,nbre_iteration,A,B,p0)
   
  nbre_base=4;
  Ai=A; Bi=B; p0i=p0;
  Aevol=zeros(nbre_etat,nbre_iteration);
  Bevol=zeros(nbre_base*nbre_etat,nbre_iteration);
  n_adn=length(adn);
  l1=zeros(nbre_etat,n_adn);
  
  for i=1:nbre_iteration,
    
    [Ai,Bi,p0i,region]=Estimation(Ai,Bi,p0i,adn,nbre_base,nbre_etat);
    
    Aevol(:,i)=diag(Ai);
    
    Bevol(:,i)=matrix(Bi',nbre_etat*nbre_base,1);
    
  end;
endfunction;

// Estimation des paramètres (NE PAS MODIFIER)
// A chaque itération, le programme principal
// fait appel à cette fonction, qui se charge
// d'appeler les fonctions de prévision, de filtrage
// et de lissage
function [Aestim,Bestim,p0estim,l1]=Estimation(A,B,p0,adn,nbre_base,nbre_etat)
  Aestim=zeros(A);
  Bestim=zeros(B);
  [prediction,filtering]=PrevFilt(A,B,p0,adn,nbre_base,nbre_etat);
  [smooth1,smooth2]=Smooth(A,prediction,filtering,adn,nbre_base,nbre_etat);
  for x=1:nbre_base,
	Bestim(:,x)=sum(smooth1(:,(adn==x)),'c');
  end;
  Bestim=Bestim./(sum(Bestim,'c')*ones(1,nbre_base));
  p0estim=sum(smooth1,'c')/sum(smooth1);
  z=sum(smooth2,'c');
  z=matrix(z,nbre_etat,nbre_etat);
  Aestim=z./(sum(z,'c')*ones(1,nbre_etat));
  l1=smooth1;
endfunction;

// NE PAS MODIFIER (Fonction résolvant les équations de prévision et de filtrage en parallèle)
function [prediction,filtering]=PrevFilt(A,B,p0,adn,nbre_base,nbre_etat);
  n=length(adn);
  filtering=zeros(nbre_etat,n);
  prediction=zeros(nbre_etat,n)
  z=Filtering(B,p0,adn(1));
  filtering(:,1)=z;
  for i=2:n,
	z=Prediction(A,z);
	prediction(:,i)=z;
	z=Filtering(B,z,adn(i));
	filtering(:,i)=z;
  end;
endfunction;

// NE PAS MODIFIER (Fonction résolvant les deux équations de lissage en parallèle)
function [smooth1,smooth2]=Smooth(A,prediction,filtering,adn,nbre_base,nbre_etat)
  n=length(adn);
  smooth1=zeros(nbre_etat,n);
  smooth2=zeros(nbre_etat^2,n);
  smooth1(:,n)=filtering(:,n);
  for i=n-1:-1:1,
	y=smooth1(:,i+1);
	filt=filtering(:,i);
	prev=prediction(:,i+1);
	z=Smooth1(y,A,prev,filt);
	smooth1(:,i)=z;
	smooth2(:,i+1)=Smooth2(y,A,prev,filt);
  end;
endfunction;

// Notations
// filt, prev ont même nature
// Ce sont les probas que Sn soit dans l'état i
// à l'issue du filtrage, de la prévision

// Équations de prédiction
// Renvoie la proba que Sn soit dans l'état i (Lemme 3)
function [y]=Prediction(A,filt)
  y = A' * filt;
endfunction;

// Équations de filtrage
// Renvoie la proba que Sn soit dans l'état i (Lemme 4)
function [y]=Filtering(B,prev,base)
  y = (B(:, base) .* prev) / sum(B(:, base) .* prev);
endfunction;

// Équations de lissage
// Deuxième équation de lissage du lemme 5
function [y1]=Smooth1(y0,A,prev,filt)
  y1 = filt .* (A * (y0 ./prev));
endfunction;

// Première équation de lissage du lemme 5
function [y2]=Smooth2(y1,A,prev,filt)
  y2 = A .* (filt * (y1 ./prev)');
  y2 = y2(:);
endfunction;

// Pour les données courtes
adn=lecture_adn("seq_lambda.txt");
A_init = [0.9 0.1; 0.1 0.9];
B_init = [0.1 0.3 0.2 0.4; 0.15 0.3 0.25 0.3];
p0_init = [0.3; 0.7];
nbre_etat = 2;
nbre_iteration = 100;

[region, p0, Aevol, Bevol] = result(adn, nbre_etat, nbre_iteration, A_init, B_init, p0_init);

// Affichier les résultats
disp(p0, 'The pobability:');
disp(Aevol(:, nbre_iteration), 'The final value of Aii');
disp(Bevol(:, nbre_iteration), 'The final value of Bij');

// Afficher les régions homogènes
clf();
homogene = zeros(1, length(adn));
homogene(find(region(1, :) >= 0.5)) = 1;
plot2d(1:length(adn), homogene, rect=[0, 0, length(adn), 1.2]);
xtitle("The homogene region");
halt();

// Afficher l'évaluation de matrice A
clf();
subplot(1, 2, 1);
plot(Aevol(1, :));
xtitle("The value of A11");
subplot(1, 2, 2);
plot(Aevol(2, :));
xtitle("The value of A22");
halt();

// Afficher l'évaluation de matrice B
clf();
for i=1:8;
  subplot(2, 4, i);
  plot(Bevol(i, :))
  xtitle(sprintf("The %dth parameter of B", i));
end;

