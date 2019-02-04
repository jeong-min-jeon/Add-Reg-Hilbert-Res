% This is a Matlab code

% Load the matlab package downloaded from https://github.com/functionaldata/PACE_matlab
addpath(genpath('C:/Downloads/release2.17/')); % Path of the package

% Define necessary objects
Time=277;
n=107;
base=(1:Time)';
T=repmat(base,[n,1]);
Ni=Time*ones(1,n);
Tcell=mat2cell(T',1,Ni);
N=108;
Yhat=[];
new_mu=[];

% To run the following code, one needs to load temperature, cloudiness and electricity data in advance
for i=1:N;
disp(i);
    % Collect predictors as X
    X=horzcat(temperature,cloudiness);
    newX=X(i,:);
    X(i,:)=[];
    % Response
    Y=electricity;
    Y(i,:)=[];
    Y = reshape(Y',[1,n*Time]);
    Ycell=mat2cell(Y,1,Ni);
    % Fit
    p=setOptions('ntest1', n, 'verbose', 'off', 'bwmu_gcv', 0);
    [S,W]=FQR(Ycell, Tcell, X, p);
    % Predict
    rhohat=S{1,4};
    mu0hat=S{1,8};
    k=size(rhohat,2);
    for j=1:k
        grid=W{j,1}{1,11}(1,:);
        value=W{j,1}{1,12}(1,:);
        new_eta=W{j,1}{1,5}*newX';
        [c index1]=min(abs(grid-new_eta));
        close1a=grid(index1);
        close1b=value(index1);
        grid(index1)=[];
        value(index1)=[];
        [d index2]=min(abs(grid-new_eta));
        close2a=grid(index2);
        close2b=value(index2);
        new_mu(j)=(close2b-close1b)/(close2a-close1a)*(new_eta-close1a)+close1b;
    end
    Yhat(i,:)=mu0hat+new_mu(1:k)*rhohat';
end
