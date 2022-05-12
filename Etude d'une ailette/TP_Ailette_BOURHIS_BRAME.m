%%
%TP Ailette : Approximation numérique par différences finies
%Victor BRAME, Tristan BOURHIS, Groupe CA
%Cas 1
T0=40;
Ta=22;
L=0.5;
R=0.0028;
Hc=10;
Lambda=40;
N=30;
Dx=L/N;
x=0:Dx:L;

w=(Hc*2)/(Lambda*R);
C1=1 ;  
C2=-(2+(w*Dx^2)) ;   
C3=1 ;  

A=zeros(N-1,N-1);
A= C1*diag(ones(1,N-2),-1)+C2*diag(ones(1,N-1),0)+C3*diag(ones(1,N-2),1);

B=zeros(N-1,1);
B=[(-w*(Dx^2)*Ta)-T0;-w*(Dx^2)*Ta*ones(N-3,1);(-w*(Dx^2)*Ta)-Ta];

T=zeros(N+1,1);
T(1,1)=T0;
T(N+1,1)=Ta;
T(2:N,1)=inv(A)*B;

plot (x,T);
%%
%solution exacte
w=sqrt((Hc*2)/(Lambda*R));

A=((-exp(-w*L))/(exp(w*L)-exp(-w*L)))*(T0-Ta);
B=(exp(w*L)/(exp(w*L)-exp(-w*L)))*(T0-Ta);

T= A*exp(w*x) + B*exp(-w*x)+Ta;

plot(x,T)
%%
%Cas 2
T0=40;
Ta=22;
L=0.5;
R=0.0028;
Hc=10;
Lambda=40;
N=30;
Dx=(L/N);
x=0:Dx:L;
w=(Hc*2)/(Lambda*R);

C1=1 ;   
C2=-(2+(w*Dx^2)) ;   
C3=1 ;   

A=zeros(N-1,N-1);
A=C1*diag(ones(1,N-2),-1)+C2*diag(ones(1,N-1),0)+C3*diag(ones(1,N-2),1);
A1=zeros(N,N);
A1=[A,[zeros(N-2,1);1];[zeros(1,N-2),-1],1];

B1=zeros(N,1);
B1=[(-w*(Dx^2)*Ta)-T0;-w*Dx^2*Ta*ones(N-2,1);0];

T1=zeros(N+1,1);
T1(1,1)=T0;
T1(2:N+1)=inv(A1)*B1;

plot(x,T1)

%%
%Solution exacte
w=sqrt((Hc*2)/(Lambda*R));

A=(exp(-w*L)/(exp(w*L)+exp(-w*L)))*(T0-Ta);
B=(exp(w*L)/(exp(w*L)+exp(-w*L)))*(T0-Ta);

T1= A*exp(w*x) + B*exp(-w*x)+Ta;
plot(x,T1)
%%
%Cas 3
T0=40;
Ta=22;
L=0.5;
R=0.0028;
Hc=10;
Lambda=40;
N=30;
Dx=(L/N);
x=0:Dx:L;
w=(Hc*2)/(Lambda*R);

C1=1 ; 
C2=-(2+(w*Dx^2)) ;  
C3=1 ;   

A=zeros(N-1,N-1);
A= C1*diag(ones(1,N-2),-1)+C2*diag(ones(1,N-1),0)+C3*diag(ones(1,N-2),1);
A2=zeros(N,N);
A2=[A,[zeros(N-2,1);1];[zeros(1,N-2),-Lambda],Dx*Hc+Lambda];


B2=zeros(N,1);
B2=[(-w*(Dx^2)*Ta)-T0;-w*Dx^2*Ta*ones(N-2,1);Dx*Hc*Ta];

T2=zeros(N+1,1);
T2(1,1)=T0;
T2(2:N+1)=inv(A2)*B2;

plot(x,T2)
%%
%Solution exacte
w=sqrt((Hc*2)/(Lambda*R));

A=((exp(-w*L)*(Lambda*w-Hc))/(Lambda*w*cosh(w*L)+Hc*sinh(w*L)))*(T0-Ta)/2;
B=((exp(w*L)*(Lambda*w+Hc))/(Lambda*w*cosh(w*L)+Hc*sinh(w*L)))*(T0-Ta)/2;

T2=A*exp(w*x) + B*exp(-w*x)+Ta;

plot(x,T2)