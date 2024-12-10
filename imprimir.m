function Imprimir(S,nV,frec1,frec2,paso)
% Datos a ingresar (Parametros S, numero de puertos, Frecuencia 1, frecuencia dos)
for i = 1:size(S,1)
    for j = 1:size(S,1)
        S_sym(i,j) = sym(sprintf('S%d%d', i, j)); %Matriz Sij
    end
end

nombre = inputdlg({'Nombre'},'Guardar como', [1 50]);
file = [nombre{1}, '.s2p'];

fileID = fopen(file, 'w');

if fileID == -1     % Comprobar si el archivo se abrió correctamente
    error('No se pudo abrir el archivo para escribir.');
end

% Escribir el encabezado del archivo S2P
fprintf(fileID, ['! Ports: ' num2str(nV) ', Freqs: ' num2str(((frec2-frec1)/paso)) '\n']);
fprintf(fileID, '# Hz S RI R 50\n'); 
fprintf(fileID, '! Freq'); 
for i = 1:length(S(:,:,1))
    for j = 1:length(S(:,:,1))
fprintf(fileID, ['                  ' char(S_sym(i,j))]);
    end
end

    fprintf(fileID, '\n'); 
for k = 0:(frec2-frec1)/paso
    fprintf(fileID, num2str(frec1+k*paso)); 
for i = 1:length(S(:,:,1))
    for j = 1:length(S(:,:,1))
fprintf(fileID, ['          ' num2str(real(S(i,j,k+1))) '          ' num2str(imag(S(i,j,k+1)))]);
    end
end
fprintf(fileID, '\n'); 
end

end
