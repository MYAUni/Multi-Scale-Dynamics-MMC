clc;        % Clear the Command Window
clear;      % Clear all variables
close all;  % Close all figure windows
%%SPVF bladder cancer model
% supervised by Dr, Ophir Nave
%% Basic paramters definition
n=5;%dimension
N=100000;
%randomize N vectors
%A biologically feasible range is considered for each variable.
randVec=[];
randVec(1,:)=randi([0 1130],1,N);%M
randVec(2,:)=randi([0 10^9],1,N);%T
randVec(3,:)=randi([0 10^5],1,N);%D
randVec(4,:)=randi([0 10^6],1,N);%E
randVec(5,:)=randi([0 10^4],1,N);%R
format long 
%format shortEng
%format compact
%% apply SPVF algorithm
subF=defineF(randVec);%substitue the random vectors in F
avgV=(1/N).*(sum(subF,2));%calculate average of random vectors
sizeS=norm(avgV);%calculate the norm of the avg vector
normF= vecnorm(subF);%calculate norm of each vector in subF 
indexF=find(normF>sizeS);%find the indices in which the norm is larger than the avg
sizeF=size(indexF,2);
K=[];%set an empty matrix K
for i=1:sizeF%subset of randV that contains vectors according to indexF
    K(:,i) = randVec(:,indexF(i));
end
cs=K;%control set
p=floor((size(cs,2))./n);
% definition of A&B
B=zeros(n,n,p);

for i=1:p%set B - basis of K  
    for h=1:n
     B(:,h,i)=cs(:,(i-1)*n+h);
    end
     B(:,n,i)=cs(:,i.*n);
end
A=B;%set a as matrices of B
%calculates the avg of the determinants of A
ADet=zeros(1,1,p);
for j = 1:p
  ADet(:,:,j) = abs(det(A(:,:,j)));
end
avgADet=(1/p).*sum(ADet);

%define a subset of A (determinants larger than the avg)
thresholdA=avgADet;
indexA=find(ADet>thresholdA);%find the indices where the determinants are bigger than the avg
sizeA=size(indexA,1);
subA=zeros(n,n,sizeA);
for l=1:sizeA
    subA(:,:,l) = A(:,:,indexA(l));
end
subB=subA;
%define T as a subset of subB substituted in  F
sizeT=size(subB,3);
T=zeros(n,n,sizeT);
t=defineF(subB);

for q=1:sizeT
    for g=1:n
     T(:,g,q)=t(:,(q-1)*n+g);
    end
     T(:,n,q)=t(:,q.*n);
end
R=[];
%conjugate transpose
ctT=zeros(n,n,sizeT);
for z=1:sizeT
    R(:,:,z)=T(:,:,z)';
    ctT(:,:,z)=T(:,:,z)*R(:,:,z);
end
%find eigenvalues and eigenvectors for each matrix in T
eT=[];
%ee=[];
for e=1:sizeT
  %  ee(:,:,e)=(eig(ctT(:,:,e)));
    eT(:,:,e)=sqrt(eig(ctT(:,:,e)));
end
sortedET=sort(eT,1);
sizeeTM=size(sortedET,3);
maxGap=zeros(1,1,sizeeTM);
for y=1:sizeeTM
    for w=2:n
        if (maxGap(:,:,y)<abs(sortedET(w,:,y)./sortedET(w-1,:,y)))
            maxGap(:,:,y)=max( maxGap(:,:,y),abs(sortedET(w,:,y)./sortedET(w-1,:,y)));
        end
    end
end
[maxVal,maxIndex] = max(maxGap(:,:,:));%find the max gap for all matrices
[maxEigVec, maxEigVal]= eig( ctT(:,:,maxIndex));%the correspinding eigenvectors


disp('The eigenvectors are:');
disp(maxEigVec);

disp('The eigenvalues are:');
disp(maxEigVal);

%For later calculation
%TsA=maxEigVec'\ eye(size(maxEigVec));
%TsAa=TsA;

%maxEigVeca=maxEigVec';

function F=defineF(x)
% Parameters of bladder cancer model
a=75;  
b=10^5;
r=0.023;
m=28740;
mu1=21.05; 
mu2=0.347;
mu3=0.2;
mu4=0.72;
p1=2.3;
p2=0.0125;
p3=(3.7)*10^(-6);
p4=1.44*10^(-5);
d0=1.032*10^3;
gamma=9.12;
beta= 5.2;
eta=7.2*10^(-2);
k=10^9;



F(1,:)=-mu1.*x(1,:)+m; %Variable M
F(2,:)=r.*(1-x(1,:)./(x(1,:)+a)).*x(2,:).*(1-x(2,:)./k)-p1.*x(2,:).*x(1,:)./(x(1,:)+a)-p3.*x(4,:).*x(2,:).*exp((-x(5,:).*x(2,:)./k).*(1-x(1,:)./(x(1,:)+a)));%Variable T
F(3,:)=d0-mu2.*x(3,:)+beta.*(1-x(3,:)./b).*x(2,:).*p1.*x(1,:)./(x(1,:)+a);%Variable D
F(4,:)=gamma.*x(3,:)-mu3.*x(4,:)-p4.*x(4,:).*x(5,:);%Variable E
F(5,:)=eta.*x(3,:).*(1-x(1,:)./(x(1,:)+a))-mu4.*x(5,:)-p2.*x(5,:).*x(1,:)./(x(1,:)+a);%Variable R
end

