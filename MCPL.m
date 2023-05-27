function F = MCPL(dataname,X,V,n,c,k,alpha,lambda,beta,gamma,miu)

%% =========================Initialization==========================
% Gv Fv
for iv = 1:V
    S = constructS_PNG(X{iv}, k, 1); %% input dv * n
    S1= (S+S')/2;
    
    Gv{iv}=S1; 
    
    D1 = diag(sum(S1));
    Lv{iv} = D1-S1;
    [Fviv, ~, ~]=eig1(Lv{iv}, c, 0);
    Fv{iv} = Fviv;
end

% Mv dv pv qv rv
for iv=1:V
    dv(iv) = size(X{iv},1); % dimension of each view
    Mv{iv} = ones(dv(iv), c);
    pv(iv) = 1.0/V;
    qv(iv) = 1.0/V;
    rv(iv) = 1.0/V;
end

% Qv LGv
for iv = 1:V
    Qv{iv} = ones(c, c);
    
    D = diag(sum(Gv{iv}));
    LGv{iv} = D-Gv{iv}; 
end

% A
A = ones(n, n);

% F H
F = ones(n, c);
F = orth(F);

H = ones(n, c);
H = orth(H);

%% =========================Normalization=========================
% Normalization
for i = 1:V
    fea = X{i}';
    fea = mapminmax(fea,0,1);  % n * dv
    X{i}= fea'; % dv * n
end

%% ==========================Optimization===========================
iter = 1;
maxIter = 5;
obj = zeros(1,maxIter);
objfile = strcat('result/',dataname,'Obj.csv');% record obj value
fp = fopen(objfile, 'a', 'n', 'utf8');
while iter<=maxIter
    
    % Updating pv,qv,rv
    for iv = 1:V
        hv(iv) = norm(X{iv}'* Mv{iv}-H,'fro');
    end
    
    for iv = 1:V
        qv(iv) = 0.5/sqrt(norm(Gv{iv}-A,'fro')^2+trace(H'*LGv{iv}*H));
        rv(iv) = 0.5/norm(Fv{iv}*Qv{iv} + H - F,'fro');
        pv(iv) = hv(iv)/sum(hv);
    end
    
    % Updating Mv
    for iv = 1:V
        I = eye(dv(iv));
        Mv{iv} = 1.0/pv(iv) * inv(1.0/pv(iv)*X{iv}*X{iv}'+alpha*I)*X{iv}*H;
    
        Mv{iv}(Mv{iv}<0) = 0;
    end
    
    % Updating A
    for i=1:n
        for j=1:n
            Q(i,j)=0.5*norm(F(i,:)-F(j,:))^2;
        end
    end
    
    up = zeros(n,n);
    down = 0;
    for iv = 1:V
        up = up + qv(iv)*Gv{iv}';
        down = down+qv(iv);
    end
    A = (-0.5*Q+beta*(H*H')+lambda*up)/(lambda*down+beta);
    
    A =(abs(A)+abs(A'))/2.0;
    
    % Updating H
    %%%L 
    In = eye(n);
    L1 = beta*(-2*A'+In);
    
    sumLGv = zeros(n, n);
    for iv = 1:V
        sumLGv = sumLGv + qv(iv) * LGv{iv};
    end
    
    L = L1+miu *sumLGv;
    
    [~,y] = eig(L);
    m = diag(y);
    lamdaMax = max(m); 
    
    %%%C
    C1 = zeros(n, c);
    for iv = 1:V
        C1 = C1 + 1.0/pv(iv)*X{iv}'*Mv{iv};
    end
    C2 = zeros(n, c);
    for iv = 1:V
        C2 = C2 + rv(iv)*(F-Fv{iv}*Qv{iv});
    end
    
    C = C1+gamma*C2;

    %%%calculate
    H1 = ones(n,c);
    I=eye(n);
    while 1
        M = 2*(lamdaMax*I-L)*H+2*C;
        [UU,~,VV] = svds(M);
        H = UU*VV';

        if(norm(H-H1))<1e-3
            break;
        end
        H1 = H;
    end
    
    % Updating F
    D = diag(sum(A));
    La = D-A; 
    L = La; 

    [~,y] = eig(L);
    m = diag(y);
    lamdaMax = max(m); 

    C1 = zeros(n, c);
    for iv = 1:V
        C1 = C1 + rv(iv)*(Fv{iv}*Qv{iv}+H);
    end
    C = gamma*C1;

    F1 = ones(n,c);
    I=eye(n);
    while 1
        M  = 2*(lamdaMax*I-L)*F+2*C;
        [UU,~,VV] = svds(M);
        F = UU*VV';
        
        if(norm(F-F1))<1e-3
            break;
        end
        F1 = F;
    end
    
    % Updating  Qv
    for iv=1:V
        Qv{iv} = inv(Fv{iv}'*Fv{iv})*(Fv{iv}'*F-Fv{iv}'*H);
    end
    
    % calculating obj value
    tempObj = 0;  % Item 1
    for iv = 1:V
        tempObj = tempObj + 1.0 / pv(iv) * norm(X{iv}' * Mv{iv} - H, 'fro')^2;
    end
    
    sumMv = 0;  % Item 2
    for iv = 1:V
        sumMv = sumMv + norm(Mv{iv},'fro')^2;
    end
    tempObj = tempObj + alpha * sumMv;
    
    sumGA = 0;  % Item 3
    for iv = 1:V
        sumGA = sumGA + qv(iv) * norm(Gv{iv} - A, 'fro')^2;
    end
    tempObj = tempObj + lambda * sumGA;
    
    tempObj = tempObj + beta * norm(A-H* H', 'fro')^2;% Item 4
    
    sumFHF = 0;  % Item 5
    for iv = 1:V
        sumFHF = sumFHF + rv(iv) * norm(Fv{iv}*Qv{iv} + H - F, 'fro')^2;
    end
    tempObj = tempObj + gamma * sumFHF;
        
    D = diag(sum(A));
    La = D-A; 
    tempObj = tempObj + trace(F'*La*F); %6
        
    sumFLGvF = 0;
    for iv = 1:V
        sumFLGvF = sumFLGvF + qv(iv)*trace(F'*LGv{iv}*F);
    end
    tempObj = tempObj + miu * sumFLGvF; %7
    
    obj(iter) = tempObj;
    fprintf(fp, '%.7f,%.7f,%.7f,%.7f,%.7f,%0.6f\n', alpha,lambda,beta,gamma,miu,obj(iter)); % 一行两个数据，用逗号分隔；每行结束后加上\n换行

    % convergence checking
    if iter>1 && abs(obj(iter)-obj(iter-1))/obj(iter-1) < 1e-3 %
        break;
    end
    iter = iter+1;
end
fprintf(fp, '\n');
fclose(fp);

%% ==========================Plot===========================

%plot(obj,'color','b');  % plot object function value according to iterations
%title(dataname);
%xlabel('Nummber of iterations');
%ylabel('Object function value');

end