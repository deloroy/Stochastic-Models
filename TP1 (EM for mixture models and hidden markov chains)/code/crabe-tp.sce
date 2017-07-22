clear

// Loi normale
function [x]=normale(y,m,s2)
  x=%e^(-(y-m).^2/2/s2)/sqrt(2*%pi*s2)
endfunction;

// Test du Chi2
function [q] = test_chi2(x, normale)
  f_normale = normale / sum(normale);
  fx = x / sum(x);
  chi2 = sum(x) * sum((f_normale - fx).^2 ./ f_normale);
  [p, q] = cdfchi("PQ", chi2, length(fx) - 1);
endfunction;

// Ouvrir le fichier de données (nombre de crabes par intervalle)
x=fscanfMat('crabe.txt');
x=x';


// intervalles
y=.580+.002+.004*[0:28];
yM=y+.002;
ym=y-.002;
Max=25;

// La fonction empirique
Fx = cumsum(x) / sum(x);
clf();
plot(y, Fx);
xgrid();
xtitle("The empiric function");
halt();

// La fonction de densité
fx = x / sum(x);
clf();
plot(y, fx);
xgrid();
xtitle("The probability density function");
halt();

// Dessiner la loi normale correspondante
x_mu = sum(x .* y) / sum(x);
x_var = sum((x.*(y - x_mu)^2)) / sum(x);
x_normale = normale(y, x_mu, x_var);
clf();
xgrid();
xtitle("The gaussian distribution");
plot2d(y, x_normale);

// Tracer l'histogramme
x_hist = x / (sum(x) * 0.004)
bar(y, x_hist);
halt();


// Test du chi2
p1 = test_chi2(x, x_normale);


// Données
pi0=[1; 3]/2/2;
pi=pi0;
mu=[.57; .67];
s2=[1; 1]/10000;

rho=ones(2, 1000);

// Algorithme EM pour les crabes avec deux populations
//------------------------------

N=1000;
R=zeros(5,N+1);
R(:,1)=[mu(1);mu(2);pi(1);s2(1);s2(2)];

z=[];
for i=1:29;
  z=[z, ones(1, x(i)) * y(i)];
end

for k=1:N;
  for i=1:2;
    rho(i,:) = pi(i) * normale(z, mu(i), s2(i)) ./ (pi(1) * normale(z, mu(1), s2(1)) + pi(2) * normale(z, mu(2), s2(2)));
  end;

  for i=1:2;
    mu(i) = sum(rho(i,:) .* z) / sum(rho(i,:));
    s2(i) = sum(rho(i,:) .* (z - mu(i))^2) / sum(rho(i,:));
    pi(i) = sum(rho(i,:)) / length(z);
  end;

  R(:, k+1) = [mu(1);mu(2);pi(1);s2(1);s2(2)];
end;

// Affichages
n1 = normale(y, mu(1), s2(1));
n2 = normale(y, mu(2), s2(2));
normale_2 = pi(1)*n1 + pi(2)*n2;
clf();
bar(y, x/(sum(x)*0.004));
plot2d(y, normale_2);
plot2d(y, pi(1)*n1,style=[color("red")]);
plot2d(y, pi(2)*n2,style=[color("green")]);
xgrid();
xtitle("The estimation with two guassian distributions");
halt();

// Tracer l'évaluation des paramètres mu
clf();
plot2d(0:length(R(1,:))-1, [R(1,:)', R(2,:)'], [1, 2]);
xtitle("The evaluation of mu");
legends(["mu_1", "mu_2"], [1, 2]);
halt();

// Tracer l'évaluation des paramètres s2
clf();
plot2d(0:length(R(4,:))-1, [R(4,:)', R(5,:)'], [1, 2]);
xtitle("The evaluation of var");
legends(["var_1", "var_2"], [1, 2]);
halt();

// Tracer l'évaluation du paramètre p1
clf();
plot2d(R(3,:));
xtitle("The evaluation of p1");
halt();

// Test du chi2
p2 = test_chi2(x, normale_2);




// Trois populations

pi0=[1; 2; 1] / 4;
pi=pi0;
mu=[0.6; 0.63; 0.67];
s2=[1; 1; 1]/10000;

rho=ones(3, 1000);

// Algorithme EM pour les crabes avec trois populations
//------------------------------

N=1000;
R=zeros(8,N+1);
R(:,1)=[mu(1);mu(2);mu(3);pi(1);pi(2);s2(1);s2(2);s2(3)];

for k=1:N;
  for i=1:3;
    rho(i,:) = pi(i) * normale(z, mu(i), s2(i)) ./ (pi(1) * normale(z, mu(1), s2(1)) + pi(2) * normale(z, mu(2), s2(2)) + pi(3) * normale(z, mu(3), s2(3)));
  end;

  for i=1:3;
    mu(i) = sum(rho(i,:) .* z) / sum(rho(i,:));
    s2(i) = sum(rho(i,:) .* (z - mu(i))^2) / sum(rho(i,:));
    pi(i) = sum(rho(i,:)) / length(z);
  end;

  R(:, k+1) = [mu(1);mu(2);mu(3);pi(1);pi(2);s2(1);s2(2);s2(3)];
end;

// Affichages
n1 = normale(y, mu(1), s2(1));
n2 = normale(y, mu(2), s2(2));
n3 = normale(y, mu(3), s2(3));
normale_3 = pi(1)*n1 + pi(2)*n2 + pi(3)*n3;
clf();
bar(y, x/(sum(x)*0.004));
plot2d(y, pi(1)*n1,style=[color("red")]);
plot2d(y, pi(2)*n2,style=[color("green")]);
plot2d(y, pi(3)*n3,style=[color("orange")]);
plot2d(y, normale_3);
xgrid();
xtitle("The estimation with three guassian distributions");
halt();

// Test du Chi2
p3 = test_chi2(x, normale_3);

// Tracer l'évaluation des paramètres mu
clf();
plot2d(0:length(R(1,:))-1, [R(1,:)', R(2,:)', R(3,:)'], [1, 2, 3]);
xtitle("The evaluation of mu");
legends(["mu_1", "mu_2", "mu_3"], [1, 2, 3]);
halt();

// Tracer l'évaluation des paramètres s2
clf();
plot2d(0:length(R(6,:))-1, [R(6,:)', R(7,:)', R(8,:)'], [1, 2, 3]);
xtitle("The evaluation of var");
legends(["var_1", "var_2", "var_3"], [1, 2, 3]);
halt();

// Tracer l'évaluation du paramètre p1
clf();
plot2d(0:length(R(4,:))-1, [R(4,:)', R(5,:)'], [1, 2]);
xtitle("The evaluation of p");
legends(["p1", "p2"], [1, 2]);
