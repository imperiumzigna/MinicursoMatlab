% exemplo fun��o 01
function [Autovalores,Autovetores] = exemplo_funcao01(M,erro)
% JACOBI
% Descri��o: Solu��o do problema de autovalores pelo m�todo de Jacobi
%            consistindo de transforma��es para triangularizar a matriz M
% Autovalores = vetor dos autovalores
% Autovetores = matriz cujas colunas sao os autovetores
% n = tamanho da matriz
% M = matriz para qual se quer determinar os autovalores e autovetores
% V = matriz unitaria para tranforma�ao
%
n = length(M);  % dimens�o da matriz A
Autovetores = eye(n); % prealoca��o da matriz unit�ria de dimens�o n (apenas para maior velocidade)

% Erro a ser considerado
if erro==0
erro = 1.0*10^(-6);
end
% determina�ao do maior coeficiente fora da diagonal prncipal de A
% ir e ic s�o a posi��o i,j do maior elemento fora da diagonal
itm = 200;  % n�mero m�ximo de itera��es
it = 0;
while it <= itm
    t=0;
    m=n-1;      
    for i=1:m
        j1=i+1;
        for j=j1:n
            if abs(M(i,j)) > t
                t=abs(M(i,j));
                ir=i;
                ic=j;
            end
        end
    end
    if it==0
       t1=t*erro; 
    end
    if t > t1
        % calculo do angulo de rota�ao
        ps = M(ir,ir) - M(ic,ic);
        % ta = tangente do angulo de rota��o
        ta = (-ps + sqrt(ps*ps+4*t*t)) / (2*M(ir,ic));
        % calculo de  cos(ta) do angulo de rota�ao
        c = 1./sqrt(1+ta*ta);
        % calculo de  sen (ta)  do angulo de rota�ao
        s = c*ta;
        % multiplicar a matrix de rota�ao por V armazenando em V
        for i=1:n
            p = Autovetores(i,ir);
            Autovetores(i,ir) = c*p+s*Autovetores(i,ic);
            Autovetores(i,ic) = c*Autovetores(i,ic)-s*p;
        end
        i=1;
        while (i-ir) ~= 0
            % Aplica�ao da transforma�ao ortogonal nas linhas de A
            p = M(i,ir);
            M(i,ir) = c*p + s*M(i,ic);
            M(i,ic) = c*M(i,ic) - s*p;
            i = i+1;
        end
        i = ir + 1;
        while (i-ic) ~= 0
            % Aplica�ao da transforma�ao ortogonal nas colunas de A
            p = M(ir,i);
            M(ir,i) = c*p + s*M(i,ic);
            M(i,ic) = c*M(i,ic) - s*p;
            i = i +1;
        end
        i=ic+1;
        while (i-n) <= 0
            p = M(ir,i);
            M(ir,i) = c*p + s*M(ic,i);
            M(ic,i) = c*M(ic,i) - s*p;
            i = i + 1;
        end
        p = M(ir,ir);
        M(ir,ir) = c*c*p + 2.*c*s*M(ir,ic) + s*s*M(ic,ic);
        M(ic,ic) = c*c*M(ic,ic) + s*s*p - 2.*c*s*M(ir,ic);
        M(ir,ic) = 0;
        it = it + 1;
    end
    it = it + 1;
end
Autovalores = diag(M);
disp('Autovalores');
disp(Autovalores);
disp('Autovetores');
disp(Autovetores);

% Arquivo com os resultados obtidos
fid=fopen('jacobi_resultados.txt','w');
fprintf(fid,'Autovalores:\t\tAutovetores:\n');
for i=1:length(Autovetores)
fprintf(fid,'%f\t\t%f\n',Autovalores(i),Autovetores(i));
end
fclose(fid);
end
