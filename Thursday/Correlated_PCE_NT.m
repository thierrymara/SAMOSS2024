if exist ('OCTAVE_VERSION', 'builtin')
    warning('off')
    pkg load statistics
    %pkg load gnuplot
    struct_levels_to_print(0)
end
close all
clear all
Here=pwd;
addpath([Here,'\Benchmarks'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nvar=3;
Nsample=2^9;
for p=1:Nvar
    Type{p}='Normal';
end
%Generate the associated LPTAU sample U(0,1)
%Uniform LPTAU Sequences
for k=1:Nsample
    x=LPTAU51(k,Nvar);
    U(k,:)=2*x-1;%U(0,1)
end

%Transformed in standard normal variable
for p=1:Nvar
    Z(:,p)=sqrt(2)*erfinv(U(:,p));
end

%Correlation Matrix
CC=zeros(Nvar,Nvar);
CC(1,2)=-0.5;CC(1,3)=0.0;CC(2,3)=0.8;
CC=CC+CC'+eye(Nvar,Nvar);
%Permutation matrix for CC
Mat=eye(Nvar-1,Nvar-1);
Mat_Perm=zeros(Nvar,Nvar);
Mat_Perm(end,1)=1;
Mat_Perm(1:end-1,2:end)=Mat;

%Generate correlated sample
U=chol(CC);
X=Z*U;

%The model
y = sum(X,2);

%Compute the Sobol' indices with correlated inputs
Perm=[1:Nvar];
Xpermuted=X;
C=CC;
for p=1:Nvar
    %Build BSPCE & Compute Sobol' indices
    [SI,PCE]=BSPCE4SAFEtoolbox(Z,y);
    fprintf('Unexplained amount of variance:  %5.4f\n', PCE.Res)
    %Collect Si,Si^ind,STi,STi^ind
    Si(Perm(1),:)=SI.Si(1,:);%With uncertainty
    STi(Perm(1),:)=SI.STi(1,:);
    Si_ind(Perm(end),:)=SI.Si(end,:);
    STi_ind(Perm(end),:)=SI.STi(end,:);
    %I circularly permute the X column
    Perm=[Perm(2:end),Perm(1)];
    Xpermuted=[Xpermuted(:,2:end),Xpermuted(:,1)];
    %I permute the correlation matrix accordingly
    C=(Mat_Perm)*(C*(Mat_Perm'));
    %Generate correlated sample
    U=chol(C);
    Z=Xpermuted*inv(U);
end
fprintf('\n')
disp('Full first-order Sobol'' indices:')
Si(:,2)'
disp('Total first-order Sobol'' indices:')
STi(:,2)'
disp('Independent total-order Sobol'' indices:')
Si_ind(:,2)'
disp('Independent total-order Sobol'' indices:')
STi_ind(:,2)'
