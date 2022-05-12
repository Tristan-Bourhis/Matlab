%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auteur : 
% BOURHIS Tristan
% 
% Groupe TP CA
%
% Le 6 avril 2021
%
% TP 3 Analyse numerique 3A :
% "R�solution num�rique de l'�quation des ondes"
%
% Encadrant : Mme. IBRAHIM
%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1.
% Voir fichier "Question 1. .pdf"

%%%%%%%%%%%%%%
% 2.
N=10; % initialisation de la borne sup�rieur de i
P=40; % initialisation de la borne sup�rieur de n
L=1;  % initialisation de la borne sup�rieur de l'intervalle de x
T=1; % initialisation de la borne sup�rieur de l'intervalle de t
deltaX =L/N; % initialisation du pas en espace
deltaT = T/P; % initialisation du pas en temps
x = 0 : deltaX : L ; % Nous prenons des x entre 0 et L avec un pas delta x
t = 0 : deltaT : T; % Nous prenons des t entre 0 et L avec un pas delta t

% On appelle la fonction "Sol_Exacte" qui va calculer la solution
% exacte pour les x, t correspondants :
% Double boucle for pour remplir le graphique :
for j=1:length(t)
    for i=1:length(x)
        Uexacte(i,j) = Sol_Exacte(x(i),t(j));
    end
end
surf(t,x,Uexacte)
title('Nappe repr�sentant la solution exacte du probl�me')
xlabel('Points temproels')
ylabel('Points spaciaux')
%%



%%%%%%%%%%%%%%%%%%%%%%
% 3. METHODE EXPLICITE
%%%%%%%%%%%%%%%%%%%%%%

% a), b), c) Voir fichier "question 3. a),b),c).pdf"


% d) a.
% M�thode explicite pour N=10 et P=40 :
N=10;
P=40;
L=1; % initialisation de la borne sup�rieur de l'intervalle de x
T=1; % initialisation de la borne sup�rieur de l'intervalle de t
c=2; % initialisation de la vitesse de propagation des ondes
deltaX =L/N; % initialisation du pas en espace
deltaT = T/P; % initialisation du pas en temps
x = 0 : deltaX : L ; % Nous prenons un x entre 0 et 1 avec un pas de h=0.01.
t = 0 : deltaT : T;
alpha = deltaT^2/deltaX^2;
U = zeros(N+1,P+1);
for i=1:length(x)
    U(i,1)=sin(pi*x(i)/L)+0.25*sin(10*pi*x(i)/L);
    U(i,2)=sin(pi*x(i)/L)+0.25*sin(10*pi*x(i)/L);
end

A = zeros(N-1,N-1);
A=(c^2)*alpha*diag(ones(N-2,1),1) + (2-2*(c^2)*alpha)*diag(ones(N-1,1))+(c^2)*alpha*diag(ones(N-2,1),-1);
B = -eye(N-1,N-1);

for n=2:P
    U(2:N,n+1)=A*U(2:N,n)+B*U(2:N,n-1);
end 
figure
surf(t,x,U)

% Pour ce cas-ci, nous avons une nappe tr�s semblable � la solution exacte.
% En effet, nous respectons la condition de stabilit� P > 2*N.

%   b.
% M�thode explicite pour N=50 et P=25 :
N=50;
P=25;
L=1;
T=1;
c=2;
deltaX =L/N;
deltaT = T/P;
x = 0 : deltaX : L ; % Nous prenons un x entre 0 et 1 avec un pas de h=0.01.
t = 0 : deltaT : T;
alpha = deltaT^2/deltaX^2;
U = zeros(N+1,P+1);
for i=1:length(x)
    U(i,1)=sin(pi*x(i)/L)+0.25*sin(10*pi*x(i)/L);
    U(i,2)=sin(pi*x(i)/L)+0.25*sin(10*pi*x(i)/L);
end

A = zeros(N-1,N-1);
A=(c^2)*alpha*diag(ones(N-2,1),1) + (2-2*(c^2)*alpha)*diag(ones(N-1,1))+(c^2)*alpha*diag(ones(N-2,1),-1);
B = -eye(N-1,N-1);

for n=2:P
    U(2:N,n+1)=A*U(2:N,n)+B*U(2:N,n-1);
end 
figure
surf(t,x,U)

% On remarque que cette figure est plate. En effet, cela vient du fait que
% nous n'avons pas respect� la condition de stabiit� qui dit que P > 2*N.


%%
% e)
% calcul de la solution exacte :
%N=50;
%P=100;
N=100;
P=75;
L=1;
T=1;
deltaX =L/N;
deltaT = T/P;
x = 0 : deltaX : L ; % Nous prenons un x entre 0 et 1 avec un pas de h=0.01.
t = 0 : deltaT : T;

% On appelle la fonction "Sol_Exacte" qui va calculer la solution
% exacte pour les x correspondants :
for j=1:length(t)
    for i=1:length(x)
        Uexacte(i,j) = Sol_Exacte(x(i),t(j));
    end
end

% calcul de la solution approch�e :
%N=50;
%P=100;
N=100;
P=75;
L=1;
T=1;
c=2;
deltaX =L/N;
deltaT = T/P;
x = 0 : deltaX : L ; % Nous prenons un x entre 0 et 1 avec un pas de h=0.01.
t = 0 : deltaT : T;
alpha = deltaT^2/deltaX^2;
U = zeros(N+1,P+1);
for i=1:length(x)
    U(i,1)=sin(pi*x(i)/L)+0.25*sin(10*pi*x(i)/L);
    U(i,2)=sin(pi*x(i)/L)+0.25*sin(10*pi*x(i)/L);
end

A = zeros(N-1,N-1);
A=(c^2)*alpha*diag(ones(N-2,1),1) + (2-2*(c^2)*alpha)*diag(ones(N-1,1))+(c^2)*alpha*diag(ones(N-2,1),-1);
B = -eye(N-1,N-1);

for n=2:P
    U(2:N,n+1)=A*U(2:N,n)+B*U(2:N,n-1);
end 



tprecis=0;
a=floor(tprecis/deltaT)+1;
U(1:N+1,a);
figure
plot(x,U(1:N+1,a),x,Uexacte(1:N+1,a),'r')


tprecis=0.5;
a=floor(tprecis/deltaT)+1;
figure
plot(x,U(1:N+1,a),x,Uexacte(1:N+1,a),'r')


tprecis=0.75;
a=floor(tprecis/deltaT)+1;
figure
plot(x,U(1:N+1,a),x,Uexacte(1:N+1,a),'r')


tprecis=1;
a=floor(tprecis/deltaT)+1;
figure
plot(x,U(1:N+1,a),x,Uexacte(1:N+1,a),'r')


% Pour conlure, lorsque nous respectons la condition de stabilit�, nous
% remarquons un tr�s faible �cart entre la solution exacte et la solution
% explicite. au contraire lorsque celle-ci n'est pa respect�e, nous avons
% des �cart tr�s grand.



% f) Voir fichier 'question 3. f).pdf'



%%
%%%%%%%%%%%%%%%%%%%%%%
% 4. METHODE IMPLICITE
%%%%%%%%%%%%%%%%%%%%%%

% a), b), c) 
% Voir fichier 'question 4. a), b), c).pdf'


% d) 
%N=25;
N=10;
L=1;
T=1;
c=2;
deltaX =L/N;
x = 0 : deltaX : L ; % Nous prenons un x entre 0 et 1 avec un pas de h=0.01.
M=[400,250,200,100];
for k=1:length(M)
    P=M(k);
    deltaT = 1/P;
    alpha = (deltaT^2)/(deltaX^2);
    
    t = 0 : deltaT : 1;
    U = zeros(N+1,P+1);
    
    
    for i=1:length(x)
        U(i,1)=sin(pi*x(i)/L)+0.25*sin(10*pi*x(i)/L);
        U(i,2)=sin(pi*x(i)/L)+0.25*sin(10*pi*x(i)/L);
    end
    
    % Remplissage des matrices :
    A = zeros(N-1,N-1);
    A=(-c^2)*alpha*diag(ones(N-2,1),1) + (1+2*(c^2)*alpha)*diag(ones(N-1,1))+(-c^2)*alpha*diag(ones(N-2,1),-1);
    B = 2*eye(N-1,N-1);
    C = -eye(N-1,N-1);

    for n=2:P
        U(2:N,n+1)= (A^(-1))*(B*U(2:N,n) + C*U(2:N,n-1));
    end 
     Uexacte = zeros(N+1,P+1);
     
for j=1:length(t)
    for i=1:length(x)
        Uexacte(i,j) = Sol_Exacte(x(i),t(j));
    end
end
    erreur(k) = max(max(abs(Uexacte-U))); % Calcul de l'erreur
end
    figure
    plot(log(M), log(erreur))
    title('Courbe de log(erreur) en fontion de log(deltaT) pour deltaT=[0.0025,0.004,0.005,0.01]');
    polyfit(log(M),log(erreur),1);

    
% L'ordre de convergence est donn� grace a la commande polyfit. L'ordre de
% convergence en temps est alors de 1. 


