%%%%%%
%%%%% File name: inv_DDCP.m 
%%%%% Computing Artifact
%%%%% Author: Amrina Ferdous
%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%
%
%Regarding the computing artifact we have considered \cite{latexTarantola} paper. Considering \\
%discrete data with continuous parameter case, authors of \cite{latexTarantola} use the following \\
%back propagation where they have used $z(w)$ with two iterations $k$-th and $(k+1)$-th but without \\
%discretizing the $w$. This back propagation is estimating the frontier
%\hat{z}_{k+1}(w)= z_0 + \int dw' \sum_i \sum_j C_{p_0}(w,w') G_{k}^{i}(w') (S^{-1})^{ij} \\
%[ d_0^j - g^j(\hat{\bfz}_k) + \int dw'' G_k^j(w'') . [\hat{z}_k(w'') - z_0(w'')] \\
%
%%
function zhat=inv_DDCP(w,x,d,sigma_d,sigma,theta,K,Gfun,ffun)
m=length(w); % Points of integration grid of Tarantola paper's algorithm
W=max(w); % Maximum of the grid
n= length(x); % Total number of data points
z= zeros(1,m) ; % Initial values for the frontier
%% asserts:
%Checking whether the total number of data is okay
%assert(n>0,'n ; total number of data should be positive integer.')
assert(n~=1,'Data is not a column vector.')
%Checking whether m is greater than one data set
assert(m~=1,'Need more than one data set.')
%Checking the lengths of x and d 
%assert(length(x)~=length(d),'The parameter on the surface,x and the data vector,d should be same dimensional column vector 1xn.')
%Checking whether the standard uncertainty of the data is okay
assert(sigma_d > 0,'sigma_d ; standard uncertainty of the data is positive floating number.')
%Checking whether the prior uncertainty is okay
assert(sigma > 0,'sigma; prior uncertainty is positive floating number.')
%Checking whether theta; a parameter of the prior covariance is okay
assert(theta > 0,'theta; a parameter of the prior covariance is positive floating number.')
%%
%%%%%%%%%%%%frontier
cov = covfun(sigma,theta,w); % Covariance matrix,C_p(w,w')=\sigma^2 exp[-\frac{1}{2} \frac{(w-w')^2}{\Delta^2}]
for k=1:K ;  %K= Number of iteration and it is changeable
    G = Gfun(x, z, w); %Jacobian
    for i=1:n ;
        for j=1:n ;
            G1 = cov.*G(j,:); % Matrix multiplication of the covariance  and the jacobian
            I3 = trapz(w,G1,2); % Trapzoidal integration of G1
            G2 = G(i,:)'.* I3; % Matrix multiplication the jacobian and the Trapzoidal integration of G1                  
            I4= trapz(w,G2);  % Trapzoidal integration of G2    
            % S matrix calculation
            S(i,j)=I4; %S matrix 
            if i==j
                S(i,j)=sigma_d+ I4;
            end
        end
    end
    %Inverse of S
    Sinv= inv(S); % Inverse of S matrix
    for j=1:n ;
        G3 = G(j,:).*z ; 
        I5(j) = trapz(w,G3); % Trapzoidal integration of G3
        f = ffun(x(j),w,z);
        % little g=integration of 'f'
        g(j)= trapz(w,f); % g= Trapzoidal integration of f w.r.t. w
        % Calculating part2(j)=[ d_0^j - g^j(\hat{\bfz}_k) + \int dw'' G_k^j(w'')\\
        %.[\hat{z}_k(w'')-z_0(w'')]
        part2(j)= d(j)-g(j)+I5(j);
    end
    % Calculating \hat{z}_{k+1}(w)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Let's consider all = \hat{z}_{k+1}(w)
    all = zeros(m,m);
    for i=1:n ;
        for j=1:n ;
            all = all + cov.*G(i,:).*Sinv(i,j)*part2(j) ;
        end
    end
    % Calculating the integration of \hat{z}_{k+1}(w)
    I6= trapz(w,all,2) ;
    z=I6' ;
    zhat(:,k)= I6 ;
end
end % end of the z_hat function
%%
%%%%%%%%%%%%%%%%
%covariance function, cov
function cov = covfun(sigma,theta,w)
    [W1, W2]= meshgrid(w);
    cov = sigma.^2 *exp((-0.5 * ((W1-W2).^2 )/theta)); %calculating the covariance
end
