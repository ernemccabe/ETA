function x=Cramer(C,A,X,N)

% USAGE: Calucaltes Cramers rule for given 
%
% INPUTS:
% C
% A
% X
% N
%
% OUTPUTS:
% CC 
% I
% 

% demo


% Start Calculating

% calculate denominator
denom=det(C);
for k=1:N
    for i=1:N
        for j=1:N
            CC(i,j)=C(i,j);
        end
    end
    for i=1:N
        CC(i,k)=A(i);
    end
    x(k)=det(CC)/denom;
end