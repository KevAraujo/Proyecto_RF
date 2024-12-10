function [Matriz_S] = embeding(S)
%[estructura_parametros_S, impedancia_caracteristica] = Read_print;

%% Calcular embeding

dut_nodos=1:3;
disp('Nodos del DUT:')
disp('puerto entrada (1), puerto salida (2), nodo comun (3)');

for i=1:1:length(S.Frequencies)
    Z(:,:,i)=(50*(S.Parameters(:,:,i)+speye(2)))/(speye(2)-S.Parameters(:,:,i));   %calcular parametros z
end

[imp,nodos] = solicitar_impedancias(S);
if max(nodos,[],'all')~=0
    if isempty(find(~nodos,1))
        mat_imp=[(Z(1,1,:)-Z(1,2,:));(Z(2,2,:)-Z(1,2,:));Z(1,2,:);imp];
        matriz_nodos=[1 (3);(3) 2;(3) 0;nodos];
    else
        mat_imp=[(Z(1,1,:)-Z(1,2,:));(Z(2,2,:)-Z(1,2,:));Z(1,2,:);imp];
        matriz_nodos=[1 (max(nodos,[],'all')+1);(max(nodos,[],'all')+1) 2;(max(nodos,[],'all')+1) 3;nodos];
    end

else 
    mat_imp=[(Z(1,1,:)-Z(1,2,:));(Z(2,2,:)-Z(1,2,:));Z(1,2,:)];
    matriz_nodos=[1 3;2 3;3 0];
end
    n_nodos = max(matriz_nodos,[],'all');

nz=length(mat_imp(:,:,1));
Vn = zeros(n_nodos,1);
z_matriz=cell(nz, 1);
matriz_Y=zeros(n_nodos,n_nodos,length(S.Frequencies));
matriz_y=matriz_Y;
V_nodos = zeros(2,2);

for i=1:1:2
    V_nodos(i,:) = str2num(input(['Nodos de V' num2str(i) ' '],'s'));
    np=max(V_nodos(i,:));
    Vn(np)=1;
end
pos_V=find(Vn == 0);

for frec=1:1:length(S.Frequencies)

    n_matriz=frec;
for i=1:1:nz
    z_matriz{i}=zeros(n_nodos);
    nmax = max(matriz_nodos(i,:));
    nmin = min(matriz_nodos(i,:));
    if nmin > 0
        z_matriz{i}(nmin,nmin)=[1]*(1/mat_imp(i,:,n_matriz));
        z_matriz{i}(nmin,nmax)=[-1]*(1/mat_imp(i,:,n_matriz));
        z_matriz{i}(nmax,nmin)=[-1]*(1/mat_imp(i,:,n_matriz));
        z_matriz{i}(nmax,nmax)=[1]*(1/mat_imp(i,:,n_matriz));
    else
        z_matriz{i}(nmax,nmax)=[1]*(1/mat_imp(i,:,n_matriz));
    end
    matriz_Y(:,:,n_matriz)=matriz_Y(:,:,n_matriz)+z_matriz{i};
end
%% Calculo de matriz

supr=0;
for i=1:1:(length(matriz_Y(:,:,n_matriz))-2)
    q=pos_V(i);
    supr(i)=q;
    for j=1:1:length(matriz_Y(:,:,n_matriz))
        for k=1:1:length(matriz_Y(:,:,n_matriz))
            matriz_y(j,k,n_matriz) = matriz_Y(j,k,n_matriz) - (matriz_Y(j,q,n_matriz)*matriz_Y(q,k,n_matriz))/(matriz_Y(q,q,n_matriz));
        end
    end
    matriz_Y(:,:,n_matriz)=matriz_y(:,:,n_matriz);
end



end
if supr > 0
    matriz_Y(supr, :,:) = [];
    matriz_Y(:, supr,:) = []; 
end


Matriz_S=y2s(matriz_Y);

%%
function [mat_imp,matriz_nodos] = solicitar_impedancias(S)

i=0;
n_matriz_total=length(S.Frequencies);
mat_imp=zeros(1)
matriz_nodos=zeros(1,2)
while true
    
    disp('¿Que desea ingresar? ')
    soli=input('R=1; L=2; C=3; exit=0;  ');

     if soli == 0 %salir
         break;
     end 

    i=i+1;

    if soli == 2
         x(i)=1i*2*pi*input(['L' num2str(i) ' = ']);
    else
        if soli == 3
            x(i)=-1i/(2*pi*input(['C' num2str(i) ' = ']));
        else
            x(i)=input(['R' num2str(i) ' = ']);
        end
    end
    matriz_nodos(i,:) = str2num(input('Nodos ','s'));
    for n_matriz=1:1:n_matriz_total

        if soli == 2
         z(i)=S.Frequencies(n_matriz)*x(i);
    else
        if soli == 3
            z(i)=x(i)/S.Frequencies(n_matriz);
        else
            z(i)=x(i);
        end
    end

    mat_imp(i,:,n_matriz)=z(i);   %matriz de impedancias (asignación de valores)

    end
    %n_matriz=0; % Reiniciar cuenta
end
end
end
