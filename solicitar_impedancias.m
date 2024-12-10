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
         x(i)=j*2*pi*input(['L' num2str(i) ' = ']);
    else
        if soli == 3
            x(i)=-j/(2*pi*input(['C' num2str(i) ' = ']));
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
